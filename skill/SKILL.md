---
name: scatacseq-fetch-process
description: |
  Fetch and process raw 10x Genomics single-cell ATAC-seq data from a public
  accession into a per-sample peak x cell .h5ad plus a durable fragment file,
  using the chromap Snakemake pipeline this skill ships with (the repository
  root, referred to below as $PIPE). Use when:
  (1) Turning a GEO/SRA accession (GSE / SRP / PRJNA / SRX / SRR / GSM) into
      fragments.tsv.gz + peak x cell matrices
  (2) Building or extending config/samples.tsv from an accession, including
      scoping ATAC libraries away from their GEX and mtDNA siblings
  (3) Launching the pipeline on SLURM and diagnosing controller hangs
  (4) Deciding ENA vs SRA download, resolving three-read geometry, or checking
      barcode whitelist orientation before spending a chromap run
  This is the UPSTREAM half (accession -> .h5ad). Downstream analysis
  (clustering, differential accessibility, motifs, footprinting) is separate.
  Handles both plain 10x scATAC and 10x Multiome ATAC, and picks between them
  from the reads rather than trusting the sample sheet.
---

# scATAC-seq Fetch & Process (accession → .h5ad)

This skill drives the Snakemake pipeline — it does **not** reimplement it.
Snakemake owns the DAG, resume-on-failure and SLURM right-sizing; this skill is
the runbook for the judgment steps around it (accession resolution, ATAC
scoping, launch, hang diagnosis, post-run verification).

**`$PIPE`** = the repository root: the directory holding `workflow/`,
`config/` and `run_slurm.sh`, i.e. the parent of the `skill/` folder that
contains **the real SKILL.md**.

> **Resolve the symlink before computing `$PIPE`.** This file is normally
> installed as a symlink at `~/.claude/skills/<name>/SKILL.md`. Taking "the
> parent of my folder" literally at that location yields `~/.claude/skills/`,
> which contains no `workflow/` — every later path would then be wrong. Derive
> it instead, which works whether the file is a symlink or a real copy:
> ```
> PIPE="$(dirname "$(dirname "$(readlink -f ~/.claude/skills/scatacseq-fetch-process/SKILL.md)")")"
> test -f "$PIPE/workflow/Snakefile"   # sanity-check before using it
> ```

> **Site activation.** Run from a working copy with a filled-in
> `config/config.yaml`, not a fresh clone (which ships `<PLACEHOLDER>` paths).
> Site paths live in `$PIPE/activate.local.sh`, which `run_slurm.sh` sources
> automatically — create it once from the template and no exports are needed:
> ```
> cd "$PIPE" && cp activate.local.sh.example activate.local.sh   # then edit
> ```
> It sets `SCATAC_VENV` (venv with snakemake + the SLURM executor plugin) and
> `SCATAC_EXTRA_PATH` (chromap, FastQC, sra-tools dirs). Both may still be
> overridden per-invocation as environment variables.

## What the pipeline does, in one paragraph

`chromap --preset atac` folds adapter trimming, paired-end mapping, barcode
correction, cell-level PCR deduplication and the Tn5 +4/−5 coordinate shift into
a single pass, writing the fragment file directly. **There is no trim rule, no
BAM, and no separate dedup step** — that is by design, not an omission. MACS3
then calls a *provisional* aggregate peak set, snapATAC2 builds the peak × cell
matrix and computes TSS enrichment, and a deterministic gate applies fixed
thresholds from `config.yaml`. The fragment file is the durable artifact; peaks
are contingent on the population called, so per-cluster recalling is downstream
work.

---

## How much of this is actually validated

Say this plainly to the user rather than implying the pipeline is broadly proven.

| | Completed runs | Distinct datasets |
|---|---|---|
| Plain 10x scATAC, GRCh38 | 5 | **2** (10x PBMC 500 and PBMC 5K) |
| Multiome ATAC, mm10 | 1 | 1 (GSE219015) |
| Multiome ATAC, GRCh38 | 0 | 0 |
| Non-PBMC human tissue | 0 | 0 |
| ATAC ↔ GEX pairing | 0 | 0 |

mm10 and multiome are each proven exactly once, on the same run. The five 10x
runs are subsamples, a reverse-complement check and a re-run of two public PBMC
datasets, so they are not five independent validations.

Consequences for how to advise:
- A first run on new tissue, chemistry or genome is a **test**. Recommend
  `keep_fastq: true` so the inputs survive if the output turns out wrong.
- The QC floors sit 4-5x below what good data produces and were set from two
  PBMC samples. They are known permissive, not known discriminating. Default to
  `qc.mode: warn`; suggest `strict` only once the user trusts the numbers for
  their tissue.
- Do not describe a metric as "good" on a first run without checking
  `mapping_rate` first — see the QC table in Phase 4 for why.

---

## Phase 1 — Build the sample sheet from an accession

```bash
# activate your Python env first
python workflow/scripts/prepare_runs.py GSE123456 -o config/samples.tsv
```

This is a **network step, run once**. It resolves GSE/SRP/PRJNA/GSM/SRX/SRR,
scoping correctly (a study fans out to all its samples; a run gives just itself).

### Five things it decides for you — check each one

**1. ATAC scoping.** A GEO series routinely carries ATAC, GEX and
mitochondrial-enrichment libraries side by side. Non-ATAC libraries are dropped
and the reason printed per run. If it drops everything, the study's
`LibraryStrategy` values are unusual — rerun with `--keep-non-atac` to see them
all before concluding the accession is wrong.

**2. Read geometry (the part that most often goes wrong).** 10x scATAC has
**three** reads: two ~50 bp genomic mates and one **16 bp cell barcode**.
Submitters name them inconsistently (`R1/R2/R3`, `R1/R2/I2`, `R1/I1/R2`), so
roles are assigned by **measured read length** via HTTP range requests, not by
filename. In `samples.tsv`, **`r2_url` must be the 16 bp barcode read** — this
is the single column worth checking by eye before launching.

**3. ENA vs SRA per run — expect SRA, and do not read that as a problem.**

ENA classifies the 16 bp cell-barcode read as a *technical* read and does not
mirror it. Surveyed across 4,000 public scATAC runs on ENA:

| files in ENA `fastq_ftp` | runs |
|---|---|
| 2 (both genomic reads only) | 3,815 |
| 0 or 1 | 185 |
| **3 (what ATAC needs)** | **0** |

So for scATAC the acquisition model is the REVERSE of the scRNA sibling: SRA is
the primary route and ENA the exception. `prepare_runs.py` still probes ENA
first (it costs one API call and a handful of studies do expose usable files via
`submitted_ftp`), but a sheet where every row says `source=sra` is the expected
outcome, not a degraded one.

Practical consequence: the SRA path carries essentially all real ATAC data, so
its failure modes matter more than the ENA ones. In particular `prefetch`
silently skips runs above its default 20 GB cap -- see Troubleshooting.

**4. Genome, from the organism SRA reports.** Not from a default. A study
spanning species is written as **one sheet per genome** and no combined file:

```
This study spans 2 genomes, so it was split into one runnable sheet per genome:
  config/samples.GRCh38.tsv   16 run(s), genome=GRCh38
  config/samples.mm10.tsv      2 run(s), genome=mm10
```

Launch each separately. A mixed sheet is refused at parse time, because merging
peak sets across species produces an object nobody should use.

This mattered on GSE219015, which is 16 human and 2 mouse samples. The `genome`
column used to be whatever `--genome` said, defaulting to GRCh38 for all 18. One
of the mouse samples was then aligned to the human genome: 97% of reads failed
to map, and the 3% that did produced a complete object with FRiP 0.70 and TSS
enrichment 21 that passed every QC check. Nothing in the output said "wrong
species".

**5. Modality — or `unknown`, which is a real answer.** `10x` vs `multiome`
selects the whitelist and the barcode offset. It is read from the protocol text
in SRA and GEO, looking for the kit name, a stated barcode length (24 bp means
ARC), or a paired RNA library. When none of those appear it writes `unknown`
rather than guessing `10x`, because a guess and a confirmed value would
otherwise look identical downstream.

You do not need to correct `unknown`. Before any alignment,
`detect_barcode_orientation.py` matches sampled reads against both whitelists in
both orientations and uses whichever wins. If the sheet names a modality and the
reads disagree, the run stops:

```
MODALITY MISMATCH: samples.tsv says modality=10x, but the reads
match multiome much better (0.9798 vs 0.0000).
```

---

## Phase 2 — Launch on SLURM

```bash
sbatch run_slurm.sh          # controller as a job
./run_slurm.sh -n            # or dry-run on the login node first
```

**Always dry-run first.** It parses `samples.tsv`, resolves the DAG and reports
the job count — and it is where a bad sheet fails in one second instead of an
hour.

**Smoke-test one sample** before committing a whole study, by targeting its
`.h5ad` so the aggregating rules are skipped:

```bash
sbatch run_slurm.sh <workdir>/h5ad/<SAMPLE>.h5ad
```

That pulls the full scientific path: `check_versions → download_fastq → fastqc →
detect_barcode_orientation → chromap → macs3_peaks → fragments_to_h5ad`.

Nothing needs the `big_memory` partition: chromap's whole point is a low memory
footprint (~10–20 GB on a human genome, against ~200 GB for STAR). Leave
`big_memory` for the scRNA pipeline.

---

## Disk: ~200 GB per sample on shared storage, and run one download at a time

Measured on GSE219015/SRR28197504 (739M spots), with the current settings:

| stage | size | where |
|---|---|---|
| `.sra` download | 24 GB | `/scratch` (node-local) |
| `fasterq-dump` output + sort temp | ~500 GB peak | `/scratch` (node-local) |
| FASTQ handed to chromap | **202 GB** | workdir (shared) |
| fragments + peaks + h5ad | ~3 GB | workdir (shared) |

Two settings decide how much of that survives:

| Setting | Effect |
|---|---|
| `keep_fastq: false` | FASTQ deleted once chromap consumes it — workdir keeps ~3 GB |
| `keep_fastq: true` | 202 GB stays until you delete it |
| `compress_fastq: false` | no gzip; saves ~5 h on a sample this size |
| `compress_fastq` unset | follows `keep_fastq` |

**Keep the FASTQs on any run whose correctness is not yet established**, and
delete them only after the output has been checked. They are only cheap to throw
away once you know you will not need them: re-fetching is a 24 GB download plus
an hour of extraction. Compression is a separate question — for a file you plan
to delete in a day, gzip costs hours and buys nothing.

**Staging is node-local.** `fasterq-dump` writes its output and sort temp to
`$SCATAC_STAGE_ROOT` (default `/scratch`, 1.9-3.6 TB per node), so the ~500 GB
peak never touches shared storage. Only the finished FASTQ is copied back. If
`/scratch` is missing or unwritable the rule falls back to the workdir and says
so.

Two things that had to be handled explicitly: `-t "$STAGE"` must be on **both**
the primary and the retry invocation, or a retry writes its sort temp into the
working directory; and the cleanup trap does not run on `scancel` (SIGKILL runs
no trap), so each download job also sweeps staging dirs whose job ID is no longer
in `squeue`.

The SLURM profile caps `download_fastq` at **one job at a time**
(`download_slots: 1`). Raising it multiplies both the node-local peak and the
shared-storage footprint.

**Why not stream.** `fasterq-dump --stdout --split-spot` would skip the
uncompressed stage entirely, but its headers carry no read index, so
demultiplexing would be positional and a spot with an unexpected read count would
silently misassign the barcode read. That trades a visible disk problem for an
invisible correctness one.

**If write access to a roomier filesystem becomes available**, the fix is a
one-line `workdir:` change in the config -- no code change.

**This is not ATAC-specific.** `scrnaseq-fetch-process` has the identical
download rule. It rarely bites there only because ENA serves most scRNA studies,
so `fasterq-dump` seldom runs -- but its own `samples.GSE219015.tsv` has 15
SRA-sourced rows and the same exposure.

## Phase 3 — Watch progress (don't confuse "running" with "hung")

```bash
squeue -u "$USER"                    # are the rule jobs actually running?
tail -f scatac_ctl.<jobid>.log       # controller log
```

The SLURM profile sets `slurm-status-command: squeue` deliberately. Snakemake's
default status poll uses `sacct`, which talks to `slurmdbd`; when that accounting
database is down (`Connection refused ... :6819`) the executor cannot observe
jobs finishing and the whole run hangs **even though every job succeeded**. If
you see a stalled controller with no queued jobs, that is the first thing to
check — this pipeline inherits the fix, but the symptom is worth recognising.

---

## Phase 4 — Verify a completed run

```bash
cat <workdir>/qc/qc_gate.tsv                       # the deterministic gate
cat <workdir>/qc/barcodes/<SAMPLE>.orientation.txt # barcode match rates
cat <workdir>/qc/read_qc.tsv                       # raw-read flags (non-fatal)
```

`qc_gate.tsv` is produced by the same code path the in-DAG gate enforces, so the
numbers cannot drift from the pass/fail decision.

Per-sample structure is checked automatically by `check_h5ad.py`, which asserts
shape, peak coordinates in `var`, non-NaN QC columns in `obs`, and that
`uns['fragments']` points at a file that exists.

### Reading the numbers

The "observed" column below is from a real reference run: 10x
`atac_pbmc_500_nextgem` (Next GEM v1.1, GRCh38), 10% downsample, 3.0M read
triples across 2 lanes. Use it to judge whether a new sample is normal, rather
than only whether it clears the floor.

| Metric | Gate floor | Observed on good data | If it is low |
|---|---|---|---|
| `mapping_rate` | > 0.50 | **0.968** | **Check the genome first.** Wrong species gives ~0.03 |
| `forward/revcomp_match_rate` | > 0.25 | **0.939** fwd / 0.000 rc | Wrong whitelist, or `r2_url` is not the barcode read |
| `frac_reads_in_peaks` (FRiP) | > 0.15 | **0.58** | Weak Tn5 signal, or peaks called on too few fragments |
| `median_tss_enrichment` | > 4 | **21.4** | Poor chromatin signal, or a GTF/reference contig mismatch |
| `duplicate_rate` | < 0.80 | **0.17** | Over-amplified / low-complexity library |
| `frac_mito` (median) | < 0.10 | **0.002** | Poor nuclei prep |
| `n_cells` | > 100 | **451** of a ~500-cell library | Cell-calling thresholds may need revisiting per tissue |

Note the floors sit 4-5x below what good data produces. That is deliberate for a
floor, but it also means they have not yet been validated against a genuinely
marginal sample - they are known to be permissive, not known to be discriminating.

**Read `mapping_rate` before anything else.** Every other metric here is computed
from the fragment file, which by construction contains only reads that mapped.
They describe the survivors and cannot see how many were lost. A GSE219015 run
against the wrong genome mapped 2.6% of its reads and passed every other check
with FRiP 0.70 and TSS enrichment 21, because those numbers came from the small
minority of reads that happened to stick.

A corollary worth internalising: after fixing that run, FRiP fell from 0.70 to
0.40 and TSS enrichment from 21 to 15.6. **The correct run looks worse.** The
inflated numbers were a selection effect. Do not treat falling QC numbers as a
regression without checking what changed upstream.

Reference runtime: the whole 17-job DAG took ~14 min wall clock on `dept_cpu`,
including a 3.5-min GRCh38 chromap index build from scratch (4.65 GB index).
Later runs reuse that index.

### Peak filtering — what it does and does not buy you

`filter_peaks` drops peaks on unplaced scaffolds/chrM and peaks overlapping the
ENCODE blacklist. Measured on two real PBMC datasets:

| | peaks called | kept | dropped-peak strength vs kept |
|---|---|---|---|
| `atac_pbmc_500_nextgem` (10% depth) | 35,037 | 34,772 | 1.02x |
| `atac_v1_pbmc_5k` (full depth) | 89,268 | 88,426 | 0.96x |

**Blacklisted peaks are not unusually strong, and filtering them does not
meaningfully change FRiP** — at either depth. Do not expect this step to improve
your QC numbers; it was tested for that and it does not.

Keep it anyway for interpretability: a peak inside a collapsed repeat or a
low-mappability region is not a regulatory element regardless of its height, and
should not appear in a differential-accessibility result.

---

## Phase 5 — Batch across studies

Three tools, all in `$PIPE`:

| Tool | Use |
|---|---|
| `run_batch.sh` | loop a manifest of accessions, one workdir each |
| `snakemake_status.sh` | DB-free status of a run: RUNNING / SUCCESS / FAILED / INCOMPLETE |
| `reconcile.py` | after a batch, classify every sample as done / rerun / broken |

```bash
PIPE=$PIPE ./snakemake_status.sh                # newest run in this repo
PIPE=$PIPE ./snakemake_status.sh scatac_ctl.<jobid>.log
```

`snakemake_status.sh` uses only `squeue` and the log files, never `sacct` or
`seff`. The SLURM accounting database is exactly what goes down on this cluster,
and a status tool that depends on the thing most likely to be broken is useless
at the moment you need it.

It also reports **RUNNING-BUT-LOOPING** (exit 5), which plain `squeue` cannot.
Snakemake's SLURM executor submits each rule as a job that re-invokes snakemake,
so a rule that fails and retries does so inside that nested session: `Error in
rule` and `Trying to restart` land in the *rule's* log, never the controller's.
From the controller the job just looks like it is still running. `download_fastq`
once retried five times over 15 hours, re-downloading 24 GB each cycle, while the
controller log showed nothing but "Job 4 has been submitted".

Still one study per working directory — set `workdir:` per config.

---

## Multiome — supported and verified

Ran end to end on GSE219015 (mouse, mm10): **8,755 cells, 96.8% of reads mapped,
FRiP 0.40, TSS enrichment 15.6**. The `allow_unverified_multiome` flag is no
longer required.

**Why it needs special handling.** In 10x Multiome the ATAC and GEX libraries
carry *different barcode sequences for the same gel bead*, corresponding by line
number across the two ARC whitelists. Without translation the two modalities get
disjoint `obs_names` and cannot be joined — and, crucially, the run still
*completes* and produces a plausible-looking object.

**What is built** (milestone M5):
- `resolve_arc_translation.py` pairs the two ARC whitelists line-for-line
  (736,320 pairs) and writes chromap's `--barcode-translate` table. Note the
  column order is `GEX<TAB>ATAC`: chromap *keys on column 2 and emits column 1*.
  Getting this backwards yields near-zero valid barcodes; an integration test
  (`tests/integration/verify_barcode_translate.py`) pins it at 120/120.
- `detect_barcode_orientation.py` probes all four modality x orientation
  combinations and emits an *oriented* whitelist and translation table.
- Multiome ATAC barcode reads are **24 bp** (8 bp spacer + 16 bp barcode), not
  16 — hence `MODALITY_OFFSET = {"10x": 0, "multiome": 8}` and
  `BARCODE_LENS = (16, 24)` in `select_reads.py`.

**The sheet does not decide this.** `--barcode-translate` is switched on when the
oriented translation table has content, and that table is written from the
modality the reads *measured* — empty for 10x, 736,320 pairs for multiome. A
mislabelled sheet therefore cannot skip translation and produce a well-formed,
wrongly-keyed object.

**What is still unverified:** the ATAC↔GEX join. Only ATAC libraries have been
processed, so no run has yet asserted that the two modalities share `obs_names`.
That assertion is the remaining check before trusting a paired analysis.

> **The bug this section exists to prevent.** `select_reads.py` originally
> hardcoded a 16 bp barcode while `detect_barcode_orientation.py` already knew
> about the 8 bp multiome offset. The two disagreed, so every attempt died in
> `select_reads` *after* a 25-minute download and ~40-minute extraction — five
> times, about 15 hours. Every test fixture used 10x geometry, including the
> ones named "multiome", so the same wrong assumption sat in both the code and
> its tests. When adding a modality, make at least one fixture carry that
> modality's *real* read geometry.

---

## Troubleshooting

*Entries below are only failures actually observed — from the test suite or from
real runs. Nothing here is invented.*

**"profile directory given (profiles/slurm), but no profile.yaml found"** —
observed on the first real `sbatch` submission. SLURM copies the submitted script
to a spool directory, so `${BASH_SOURCE[0]}` inside `run_slurm.sh` does NOT point
at the repo. `run_slurm.sh` now resolves the repo root as `$SCATAC_REPO` →
`$SLURM_SUBMIT_DIR` → script directory, and fails with a clear message if none of
them contains `workflow/Snakefile`. If you hit this anyway, either `sbatch` from
the repo root or `export SCATAC_REPO=/path/to/repo`.

Note this class of bug is invisible to `snakemake -n` and to running snakemake
directly — only a real `sbatch` exercises that code path.

**A run mysteriously re-does work you thought was finished, or two runs fight
over one workdir.** Observed during development, twice, both self-inflicted:

  1. *Never `rm` outputs from a workdir that has a live run against it.* Snakemake's
     lock is held by the first invocation; a second `snakemake` started elsewhere
     is a separate process and will happily write the same files. Deleting a file
     the running controller expects sends it into a retry loop while the second
     process rebuilds the same targets. If you need to force a rerun, stop the
     controller first, or use `snakemake --forcerun <rule>`.
  2. *Editing a rule's `params:` or `shell:` invalidates it and everything
     downstream*, even when the change is cosmetic for your samples. Adding an
     optional flag that renders empty for 10x still changes the rule's
     provenance, and Snakemake will re-align. Check `snakemake -n` before
     assuming an edit is free -- on a full-depth sample that is a 25-minute
     surprise.

**A metric looks unchanged after a rerun that "succeeded".** Confirm the output
file's mtime actually moved. A failed run leaves the previous file in place, and
reading it back looks exactly like a passing regression check. This produced a
false "PASS" during development.

**"Unable to guess SLURM account. Trying to proceed without."** — expected and
harmless. The profile sets `slurm-no-account: true` deliberately, because the
account probe goes through `sacct`/slurmdbd and fails exactly when that database
is down — the same outage that makes the default status polling hang.

**"MODALITY MISMATCH: samples.tsv says modality=X, but the reads match Y"** — the
sheet is wrong and the reads are right. Set `modality=Y` for that sample, or
`unknown` to let the measurement decide, and rerun. Nothing downstream ran, so
only the download is lost.

**"samples.tsv mixes genomes: ['GRCh38', 'mm10']"** — the study spans species.
One run handles one genome. Rebuild the sheet with `prepare_runs.py`, which now
writes one sheet per genome, and launch each separately.

**Mapping rate far below 0.50** — check the species before anything else.
`grep -i scientificname` on the SRA metadata for one run, and compare it to the
`genome` column. Aligning to the wrong species maps ~3% of reads and still
produces a complete, healthy-looking object.

**"Neither orientation matches the whitelist"** — the run stopped in seconds
rather than producing an empty fragment file. In order of likelihood: wrong
whitelist for the chemistry (10x scATAC needs `737K-cratac-v1.txt`;
`737K-arc-v1.txt` is Multiome only); `r2_url` is not the 16 bp barcode read; the
data is not 10x droplet scATAC.

**"chromap produced no fragments"** — check
`qc/barcodes/<sample>.orientation.txt` first, then confirm the genomic reads were
not passed in the barcode slot.

**"X has zero non-zero entries"** (from `check_h5ad`) — almost always a contig
naming mismatch between the fragments and the peaks/reference (`chr1` vs `1`).
The pipeline takes chromosome sizes from the alignment FASTA itself to prevent
this, so it points at a peak file called against a different reference.

---

## Decision: skill vs. rewriting the pipeline

Fixes belong **inside `$PIPE`**, in version control. Thresholds live in
`config.yaml` so a QC decision is reproducible and reviewable; do not adjust a
gate by hand-editing outputs or by passing one-off flags.

## Key files in `$PIPE`

| Path | What it is |
|---|---|
| `workflow/Snakefile` | the DAG |
| `config/config.example.yaml` | every tunable, with the reasoning |
| `workflow/scripts/prepare_runs.py` | accession → `samples.tsv`, ATAC scoping |
| `workflow/scripts/detect_barcode_orientation.py` | the cheap early failure |
| `workflow/scripts/fragments_to_h5ad.py` | the object everything else consumes |
| `profiles/slurm/config.yaml` | partitions, memory, squeue polling |
| `tests/` | `pytest tests/` — no cluster or network needed |

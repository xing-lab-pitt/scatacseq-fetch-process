# scatacseq-fetch-process

Turn a public **GEO or SRA accession** into per-sample **peak × cell `.h5ad`**
files plus a durable, tabix-indexed **fragment file**, on SLURM, with a
deterministic QC gate.

(Not ENA — it does not mirror the barcode read for scATAC. See below.)

Companion to [`scrnaseq-fetch-process`](https://github.com/xing-lab-pitt/scrnaseq-fetch-process),
mirroring its folder layout, multi-sample acquisition model and "pipeline
decides / skill explains" philosophy — swapping only
the three slots that differ between assays: **aligner, matrix builder, QC
metrics**.

```
accession ──> samples.tsv ──> FASTQ ──> chromap ──> fragments.tsv.gz
                                             ├──> MACS3 peaks
                                             └──> snapATAC2 matrix ──> .h5ad
```

## Status — what has actually been run

Read this before trusting a number it produces.

| | Completed runs | Distinct datasets | Notes |
|---|---|---|---|
| Plain 10x scATAC, GRCh38 | 5 | **2** | 10x PBMC 500 and PBMC 5K; the rest are subsamples, a reverse-complement check and a re-run |
| Multiome ATAC, mm10 | 1 | 1 | GSE219015: 8,755 cells, 96.8% mapped, FRiP 0.40, TSSe 15.6 |
| Multiome ATAC, GRCh38 | 0 | 0 | — |
| Non-PBMC human tissue | 0 | 0 | the 16 human samples in GSE219015 have not been run |
| ATAC ↔ GEX pairing | 0 | 0 | only ATAC libraries processed; `obs_names` never compared across modalities |

Five runs is not five independent validations — they cover two public PBMC
datasets. **mm10 and multiome are each proven exactly once, on the same run.**

Treat a first run on new tissue, chemistry or genome as a test, not a result, and
keep `keep_fastq: true` until you have looked at the output.

**The QC thresholds are permissive floors, not tuned cutoffs.** They sit 4-5×
below what good data produces and were set from two PBMC samples. They are known
to be loose; they are *not* known to catch a marginally bad sample. `qc.mode:
warn` reports without stopping — use `strict` only once you trust the numbers for
your tissue.

**One failure mode is worth knowing up front.** Every QC metric except
`mapping_rate` is computed from the fragment file, which only contains reads that
mapped. A run against the wrong genome mapped 2.6% of its reads and still
produced a complete object with FRiP 0.70 and TSS enrichment 21 — better-looking
than the correct run, because the few reads that mapped were the cleanest ones.
Check `mapping_rate` first, always.

## Design decisions worth knowing before you read the code

**There is no adapter-trim rule, no BAM, and no separate dedup step.**
`chromap --preset atac` folds trimming, mapping, barcode correction, cell-level
PCR deduplication and the Tn5 +4/−5 shift into one pass, emitting the fragment
file directly. This was checked against two reference pipelines: MAESTRO and
scATACpipe both carry a trim step **only** on their non-chromap (bwa/bowtie2)
paths — in scATACpipe, `CUTADAPT` feeds `BWA_MAP` and never chromap.
FastQC still runs, as a pure diagnostic that does not modify reads.

### Data sources: GEO + SRA only. ENA does not work for scATAC.

**GEO** for metadata, **SRA** for reads. ENA is probed but should be expected to
fail, for a specific and unavoidable reason.

| Archive | What it provides | Usable here |
|---|---|---|
| **GEO** | Series/sample metadata: which libraries exist, and which are ATAC vs GEX vs mtDNA | **Yes** — the entry point. Hosts no reads. |
| **SRA** | Raw reads as `.sra`, expanded locally with `fasterq-dump` | **Yes** — the source of essentially all data |
| **ENA** | Mirrors SRA *and* publishes ready-made FASTQ | **No** — see below |

ENA classifies the 16 bp cell-barcode read as a **technical** read and does not
publish it. Surveyed across 4,000 public scATAC runs on ENA:

| files in ENA `fastq_ftp` | runs |
|---|---|
| 2 (both genomic mates, no barcode) | 3,815 |
| 0 or 1 | 185 |
| **3 (what scATAC needs)** | **0** |

Two genomic reads with no barcode is not scATAC data — there is no way to tell
which cell anything came from. So the ENA-first model the scRNA sibling uses is
simply inapplicable here: for scRNA the barcode and cDNA are both *biological*
reads and ENA serves them; for scATAC the barcode is *technical* and ENA drops it.

Consequences worth knowing before you run anything:

- **A sheet where every row says `source=sra` is correct output**, not a
  degraded fallback. `prepare_runs.py` says so in its own output.
- **Downloads are heavier.** SRA means fetching a `.sra` and expanding it, rather
  than streaming ready-made FASTQ — roughly 3x the disk churn. Budget for it.
- **The SRA path's failure modes are the ones that matter**, because it is the
  only road. Chief among them: `prefetch` silently skips runs above its default
  20 GB cap, and most real scATAC runs exceed that. The pipeline passes an
  explicit `--max-size` and verifies a `.sra` actually appeared.
- **Integrity checking differs by source.** The `r*_md5` columns hold ENA
  checksums and are therefore almost always empty for scATAC. SRA downloads are
  covered instead by `prefetch --verify`, and by the guard that a `.sra` file
  exists before `fasterq-dump` runs. Do not read empty md5 columns as a missing
  safeguard — they are the wrong safeguard for this route.
- **A minority escape hatch:** about 4% of scATAC runs expose 3-4 files under
  ENA's `submitted_ftp` field (the submitter's original uploads rather than
  ENA's conversions). The pipeline does not read that field. It is a possible
  future optimisation for a few studies, not a general solution.

**Barcode whitelist orientation is detected, not configured.** Whether a 10x
whitelist matches the reads as written or reverse-complemented depends on the
*sequencer*, not the assay. Guessing wrong does not error — chromap completes and
emits an almost-empty fragment file. So we sample reads, count exact matches in
both orientations, and abort in seconds if neither clears a threshold. (Technique
borrowed from ENCODE_scatac's `barcode_revcomp_detect.py`.)

**Read roles are assigned by measured length, not filename.** 10x scATAC has
three reads — two genomic mates and a **16 bp** barcode — deposited under
inconsistent names (`R1/R2/R3`, `R1/R2/I2`, `R1/I1/R2`). The 16 bp read is the
barcode; that is the only stable signal.

**The peak set is provisional.** Peak calling is population-dependent, so the
aggregate MACS3 call exists to give the QC gate a FRiP number and build a first
object. Per-cluster pseudobulk recalling belongs downstream. The **fragment file
is the durable substrate** — peaks are contingent, fragments are not.

**QC thresholds live in version control**, in `config.yaml`, never chosen at
runtime. Reporting (`qc_gate`) is a separate rule from enforcement (`qc_check`)
so a failing gate never deletes its own report.

## Install

Python packages go into an existing venv with `uv` — no new conda environment:

```bash
uv pip install --python "$SCATAC_VENV" macs3 snapatac2 pysam
```

`chromap` is C++ and **not** on PyPI, so it is built from source and put on PATH:

```bash
cd <your software dir> && git clone https://github.com/haowenz/chromap.git
cd chromap && make
export SCATAC_EXTRA_PATH=<your software dir>/chromap:/opt/FastQC:<sra-tools bin>
```

The 10x scATAC whitelist `737K-cratac-v1.txt` is not bundled here; it ships with
scATACpipe (`assets/whitelist_barcodes/`) and with Cell Ranger ATAC.

`bedtools`, `bgzip` and `tabix` are deliberately **not** dependencies: fragment
compression and indexing go through `pysam`, and fragments∩peaks through
`snapATAC2`/`pyranges`.

## Quickstart

```bash
cp activate.local.sh.example activate.local.sh     # once: your venv + tool paths
cp config/config.example.yaml config/config.yaml   # then edit the paths
python workflow/scripts/prepare_runs.py GSE123456 -o config/samples.tsv
./run_slurm.sh -n                                  # dry run: check the DAG
sbatch run_slurm.sh                                # launch
```

`activate.local.sh` is git-ignored and sourced by `run_slurm.sh`. Without it the
launcher depends on whatever environment you happened to run `sbatch` from, and
a submission from the wrong shell dies two seconds later with
`MISSING on PATH: snakemake`.

`prepare_runs.py` takes a **GEO** accession (or SRP/PRJNA/GSM/SRX/SRR). Reads come
from **SRA**; expect every row of the sheet to say `source=sra` — that is normal
for scATAC, not a fallback. Set `sra_tools_bin` in the config, or put
`prefetch`/`fasterq-dump` on PATH.

It fills `genome` from the organism SRA reports, so a mouse study is never
recorded as GRCh38. **A study spanning species produces one sheet per genome**
(`samples.GRCh38.tsv`, `samples.mm10.tsv`) and no combined file, because one run
handles one genome and a mixed sheet is refused at parse time. GSE219015 is 16
human plus 2 mouse samples and splits accordingly.

### Smoke-test one sample first

Target a single sample's `.h5ad` to pull the full scientific path without the
aggregating rules:

```bash
sbatch run_slurm.sh <workdir>/h5ad/<SAMPLE>.h5ad
```

## Output

```
<workdir>/
  fragments/<sample>.fragments.tsv.gz(.tbi)   # the durable artifact
  peaks/<sample>_peaks.narrowPeak             # provisional
  h5ad/<sample>.h5ad                          # peak x cell + per-barcode QC
  combined.h5ad
  qc/qc_gate.tsv, qc/read_qc.tsv, qc/multiqc_report.html
  qc/barcodes/<sample>.orientation.txt
```

`.h5ad` contents: `X` = peak × cell counts; `var` = chrom/start/end;
`obs` = `n_fragments`, `frip`, `tss_enrichment`, `nucleosome_signal`,
`frac_mito`, `frac_dup`, `is_cell`; `uns['fragments']` = path to the fragment
file.

`uns['mapping_rate']` is the one metric not derived from the fragment file. The
fragment file only holds reads that mapped, so every other number describes the
survivors and cannot see how many were lost. Gate it with `qc.min_mapping_rate`.
A run of GSE219015 against the wrong genome mapped 2.6% of its reads and passed
every other check with FRiP 0.70 and TSS enrichment 21.

### Keeping the raw FASTQs

`keep_fastq` and `compress_fastq` are independent:

| | Effect |
|---|---|
| `keep_fastq: false` | FASTQs deleted once chromap has consumed them |
| `keep_fastq: true` | kept, so results can be checked before the inputs go |
| `compress_fastq: false` | no gzip — saves hours on a large study |
| `compress_fastq` unset | follows `keep_fastq` |

Keep them on any run whose correctness is not yet established. Uncompressed
FASTQs for one deeply-sequenced sample run to ~200 GB, but re-fetching them costs
a 24 GB download plus an hour of extraction.

## Tests

```bash
pytest tests/          # 116 tests, no cluster, no network, no real data
bash tests/test_site_config.sh   # launcher wiring; skips without activate.local.sh
```

Most of these exist because something went wrong once:

| Area | What it pins |
|---|---|
| Read roles | assignment by measured length, including the 8 bp sample-index trap and the 24 bp multiome barcode |
| Genome | follows the organism SRA reports; an explicit `--genome` that contradicts it is refused |
| Modality | `unknown` when the metadata is silent, never a guessed `10x` |
| Barcode orientation | both directions plus the abort path |
| QC gate | boundary values, missing metrics, and mapping rate |
| FASTQ flags | `keep_fastq` and `compress_fastq` stay independent |
| Sheet splitting | one runnable sheet per genome, no rows lost |
| `.h5ad` | structural checks |

Each was checked by breaking the code and confirming the tests fail. A test that
still passes with the fix reverted is not evidence.

## Multiome

Supported. In 10x Multiome the ATAC and GEX libraries carry different barcode
sequences for the same gel bead, corresponding by line number across the two ARC
whitelists. `resolve_arc_translation.py` pairs them into the table chromap needs,
and the ATAC barcode is read from bases 8-23 of a 24 bp barcode read rather than
0-15 of a 16 bp one.

Verified end to end on GSE219015 (mouse, mm10): 8,755 cells, 96.8% of reads
mapped, FRiP 0.40, TSS enrichment 15.6.

The `modality` column is a starting value, not the decision. Before any
alignment, `detect_barcode_orientation.py` matches sampled reads against both
whitelists in both orientations and stops the run if the sheet disagrees with
what it measures:

```
MODALITY MISMATCH: samples.tsv says modality=10x, but the reads
match multiome much better (0.9798 vs 0.0000).
```

Write `unknown` in that column and the measurement is used with no complaint.
`prepare_runs.py` writes `unknown` itself whenever the metadata does not say,
rather than defaulting to `10x` — a guess and a confirmed value should not look
identical downstream.

## Prior art

Each reference contributed one clearly-scoped idea:

- **scATAC-pro** (Yu et al., *Genome Biology* 2020) — the architecture template:
  an explicit split between data processing (FASTQ→matrix+QC) and downstream
  analysis. That boundary is this pipeline's DAG-vs-skill seam. Also the
  transparent FRiP + unique-fragment cell-calling philosophy.
- **MAESTRO** (Wang et al., *Genome Biology* 2020) — validates chromap as the
  default aligner and the fastq→fragments→QC→matrix flow with no trim step.
- **ENCODE_scatac** (Kundaje lab) — the `Modality` key, and the empirical
  reverse-complement detection. Note it does **not** use chromap; its own path is
  bwa/bowtie2 → BAM → dedup → fragments, which is why it needs the fastp trim
  step this pipeline does not.
- **scATACpipe** (Hu et al.) — the argument that peak calling is
  population-dependent and belongs in per-cluster pseudobulk downstream, making
  the DAG's peak set provisional.

## License

See `LICENSE`.

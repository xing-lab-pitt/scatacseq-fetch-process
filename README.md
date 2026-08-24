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

`bedtools`, `bgzip` and `tabix` are deliberately **not** dependencies: fragment
compression and indexing go through `pysam`, and fragments∩peaks through
`snapATAC2`/`pyranges`.

## Reference data — what you need and where to get it

None of this is in the repo: the files are large, and the 10x ones are covered by
10x's licence. Put them wherever you like and point `config.yaml` at them.

Budget roughly **20 GB per genome**, most of it the chromap index.

### 1. Genome FASTA + GTF

Only these two files are used from any reference bundle.

| Genome | Source used here | Notes |
|---|---|---|
| GRCh38 | 10x Cell Ranger reference **`refdata-gex-GRCh38-2024-A`** | browser download + licence acceptance; the ARC reference works equally well |
| mm10 | UCSC `mm10.fa` + `mm10.ncbiRefSeq.gtf` | direct download, below |

```bash
# mm10 from UCSC (both links verified)
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz
gunzip mm10.fa.gz mm10.ncbiRefSeq.gtf.gz     # ~2.6 GB and ~490 MB unpacked
```

**Alternatively, let a reference pipeline fetch the genome.**
[scATACpipe](https://github.com/hukai916/scATACpipe) downloads and prepares a
UCSC or Ensembl genome for you:

```
--ref_fasta_ucsc mm10            # or
--ref_fasta_ensembl homo_sapiens
```

Only the FASTA and GTF it produces are needed here — the rest of its output is
for its own workflow.

**The FASTA and the GTF must use the same contig names.** A `chr`-prefixed FASTA
against an unprefixed GTF gives zero TSS enrichment and no error message.

This is subtler than "UCSC vs Ensembl". The 10x GRCh38 reference above is *built
from* Ensembl (`Homo_sapiens.GRCh38.dna.primary_assembly.fa` + GENCODE v44) but
10x re-prefixes the contigs, so its FASTA starts `>chr1`. Dropping a raw Ensembl
FASTA next to the 10x GTF therefore breaks, even though both are "GRCh38". Check
with:

```bash
head -1 <fasta>                        # >chr1 ... or >1 ...
awk '!/^#/{print $1; exit}' <gtf>      # chr1   or 1
```

### 2. Chromap index — built, not downloaded

`resolve_chromap_index` builds it from your FASTA on first use (~4 min for
GRCh38) and records the chromap version alongside it, so an aligner upgrade
rebuilds rather than silently reusing a stale index. Expect **~12 GB per genome**.

### 3. Barcode whitelists

Which one you need depends on the chemistry, and the pipeline picks between them
from the reads:

| File | Barcodes | Used for |
|---|---|---|
| `737K-cratac-v1.txt` | 737,280 × 16 bp | plain 10x scATAC |
| `737K-arc-v1.ATAC.txt` | 736,320 × 16 bp | Multiome — the **ATAC** half |
| `737K-arc-v1.txt` | 736,320 × 16 bp | Multiome — the **GEX** half |

They ship inside the 10x software bundles: `737K-cratac-v1` with **Cell Ranger
ATAC**, and the two `737K-arc-v1` files with **Cell Ranger ARC**. Both are free
downloads from the same 10x support site.

[scATACpipe](https://github.com/hukai916/scATACpipe) mirrors the scATAC one under
`assets/whitelist_barcodes/`, as both `737K-cratac-v1.txt.gz` and a
reverse-complemented `737K-cratac-V1-rc.txt.gz`. **You only need the forward
file.** Some chemistries write the barcode reverse-complemented, which is why
that second file exists elsewhere; here `detect_barcode_orientation.py` measures
which orientation the reads use and reverse-complements the whitelist itself.
Supplying a pre-RC'd list would double-flip it.

**The two ARC files are a pair and their line order is the data.** Line *n* of the
ATAC list and line *n* of the GEX list are the same gel bead. That correspondence
is what `resolve_arc_translation.py` turns into chromap's `--barcode-translate`
table — so do not sort, deduplicate or otherwise reorder either file.

```yaml
barcodes:
  whitelist_10x:           "/refs/10x_whitelists/737K-cratac-v1.txt"
  whitelist_multiome_atac: "/refs/10x_whitelists/737K-arc-v1.ATAC.txt"
  whitelist_multiome_gex:  "/refs/10x_whitelists/737K-arc-v1.txt"
```

### 4. Blacklists

ENCODE blacklist **v2** (Amemiya, Kundaje & Boyle, *Sci Rep* 2019), from the
Boyle Lab repo. Small files, both links verified:

```bash
wget https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz
wget https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz mm10-blacklist.v2.bed.gz
```

**Check the region count, not the filename.** Several different files get called
"the mm10 blacklist", they are all real ENCODE products, and the names do not
distinguish them:

| File | Regions | Covered |
|---|---|---|
| `hg38-blacklist.v2.bed` (v2) | 636 | 227 Mb |
| `mm10-blacklist.v2.bed` (v2) | **3,435** | **239 Mb** |
| ENCFF547MET — what ENCODE_scatac uses for mm10 | 164 | 0.1 Mb |

```bash
awk '{s+=$3-$2} END{printf "%d regions, %.1f Mb\n", NR, s/1e6}' <blacklist>
```

A 164-region mm10 file is effectively no filtering. It was on this machine under
the name `mm10-blacklist.v2.bed`, which is how it went unnoticed. The v2 files
also carry an annotation column (`High Signal Region` / `Low Mappability`); the
older list is bare 3-column BED.

Optional — omit `blacklist:` from a `references:` block and filtering is skipped.
Worth knowing what it does and does not buy you: measured at two sequencing
depths, the peaks it removes were **not** disproportionately strong (1.02× and
0.96× the median). It is there for interpretability, not to improve your metrics.

### Where the pipeline expects them

Paths are per-genome in `config.yaml`, so one config can serve several genomes:

```yaml
references:
  GRCh38:
    fasta:         "/refs/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
    gtf:           "/refs/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
    chromap_index: "/refs/GRCh38.chromap"      # created on first run
    gsize: hs
    blacklist:     "/refs/blacklists/hg38-blacklist.v2.bed"
  mm10:
    fasta:         "/refs/mm10/mm10.fa"
    gtf:           "/refs/mm10/mm10.ncbiRefSeq.gtf"
    chromap_index: "/refs/mm10.chromap"
    gsize: mm
    blacklist:     "/refs/blacklists/mm10-blacklist.v2.bed"
```

`gsize` is the MACS3 effective genome size flag: `hs` for human, `mm` for mouse.
The sheet's `genome` column selects the block, and `prepare_runs.py` fills that
column from the organism SRA reports.

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

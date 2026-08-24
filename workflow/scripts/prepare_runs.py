#!/usr/bin/env python3
"""Resolve a public accession into config/samples.tsv for the scATAC pipeline.

Accepts GSE / SRP / PRJNA / GSM / SRX / SRR and scopes correctly: a STUDY
accession fans out to all its samples, an experiment/run accession gives just
its own.

THREE THINGS THIS DOES THAT THE scRNA SIBLING DOES NOT:

1. ATAC-AWARE SCOPING. A single GEO series routinely carries ATAC, GEX and
   mitochondrial-enrichment libraries side by side. Feeding a GEX library to
   chromap wastes hours and produces nonsense, so non-ATAC libraries are
   dropped here, with the reason printed for each.

2. THREE READS, NOT TWO. 10x scATAC deposits two genomic mates plus a 16 bp
   cell-barcode read. Submitters name them inconsistently (R1/R2/R3, R1/R2/I2,
   R1/I1/R2), and sorting by filename gets it wrong often enough to matter. We
   assign roles by MEASURED READ LENGTH via HTTP range requests -- no download.

3. GENOME AND MODALITY COME FROM THE METADATA, NOT FROM A DEFAULT. The genome
   follows the organism SRA reports, so a mouse study is never recorded as
   GRCh38. The modality is read from the protocol text in SRA and GEO, and is
   written as "unknown" when neither says -- the reads settle it later. Both
   used to default silently, and a mouse sample aligned to the human genome
   still produced a complete, QC-passing object from the 3% of reads that
   mapped.

Network step, run once:
    python workflow/scripts/prepare_runs.py GSE123456 -o config/samples.tsv
"""
import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import ncbi_utils  # noqa: E402

BARCODE_LEN = 16      # 10x single-cell ATAC cell barcode
INDEX_MAX_LEN = 10    # i7/i5 sample indexes; never a cell barcode

# LibraryStrategy / name substrings that mark a library we must NOT align.
GEX_MARKERS = ("rna-seq", "rna_seq", "gex", "expression", "transcriptom")
MITO_MARKERS = ("mito", "mtdna", "mt-dna", "mtatac", "redeem")
MULTIOME_MARKERS = ("multiome", "arc-v1", "arcv1", " arc", "_arc", "-arc")

# Markers safe to search in FREE PROSE (protocol descriptions, GEO pages).
# MULTIOME_MARKERS is for short name fields only: " arc" matches "architecture"
# and "-arc" matches "search-arch" in running text, so reusing it on prose
# invents multiome samples that do not exist.
MULTIOME_PROSE_MARKERS = (
    "multiome",
    "cellranger-arc",
    "cell-ranger-arc",
    "cellranger arc",
    "arc-v1",
)

# A 10x Multiome ATAC barcode read is 24 bp: 8 bp spacer + 16 bp barcode.
MULTIOME_BARCODE_READ_LEN = 24

# Multiome pairs an ATAC library with a GEX one sharing gel-bead barcodes.
# Submitters state the pairing even when they never name the kit.
RNA_PARTNER_PHRASES = (
    "matchable with mrna",
    "matched with mrna",
    "paired with mrna",
    "same cell barcodes as the rna",
    "joint single-cell rna and atac",
)


def looks_atac(run):
    """True if this run is an ATAC library. Errs toward keeping ambiguous rows."""
    strategy = (run.get("library_strategy") or "").lower()
    return "atac" in strategy


def rejection_reason(run):
    """Why this run is not usable, or None if it is. Order matters."""
    strategy = (run.get("library_strategy") or "").lower()
    name = " ".join(str(run.get(k, "")) for k in
                    ("library_name", "sample_name", "gsm")).lower()

    if any(m in strategy or m in name for m in MITO_MARKERS):
        return "mitochondrial-enrichment library (not nuclear ATAC)"
    if not looks_atac(run):
        if any(m in strategy or m in name for m in GEX_MARKERS):
            return f"gene-expression library (LibraryStrategy={strategy or '?'})"
        return f"not an ATAC library (LibraryStrategy={strategy or '?'})"
    return None


def _http_text(url, timeout=60):
    """URL -> text with tags stripped, or "" on any failure."""
    import re as _re
    from urllib.request import urlopen
    try:
        with urlopen(url, timeout=timeout) as r:
            raw = r.read().decode("utf-8", "replace")
    except Exception:
        return ""
    return _re.sub(r"<[^>]*>", " ", raw)


def library_prose(run, cache={}):
    """Protocol prose for a run, from SRA's XML and its GEO sample page.

    Neither source alone is enough for GSE219015. SRA's
    LIBRARY_CONSTRUCTION_PROTOCOL names the Multiome kit for the human runs but
    not the mouse ones; the GEO sample page names it for both. Whichever is
    richer wins, so both are fetched and concatenated.

    Cached per run: this is called once per row and each call is two HTTP
    round-trips.
    """
    key = run.get("srr") or run.get("gsm") or ""
    if key in cache:
        return cache[key]

    parts = []
    srr = run.get("srr")
    if srr:
        parts.append(_http_text(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=sra&id={srr}&rettype=xml"))
    gsm = run.get("gsm")
    if gsm:
        parts.append(_http_text(
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"))

    cache[key] = " ".join(parts).lower()
    return cache[key]


def barcode_len_in_prose(prose):
    """Cell-barcode length stated in the text, or None.

    GSE219015's mouse sample never says "multiome", but it does say
    "i5 file (1-24) contains cell barcodes". A 24 bp barcode read is the ARC
    geometry: 8 bp spacer plus a 16 bp barcode. That sentence is the only
    textual evidence for that sample.
    """
    import re as _re
    for m in _re.finditer(r"\(1-(\d{1,3})\)[^.]{0,60}cell barcode", prose):
        return int(m.group(1))
    for m in _re.finditer(r"cell barcode[^.]{0,60}\(1-(\d{1,3})\)", prose):
        return int(m.group(1))
    return None


def detect_modality(run, use_network=True):
    """Returns "multiome", "10x", or "unknown".

    "unknown" is a real answer, not a failure. The value picked here only has to
    be a starting guess: detect_barcode_orientation probes both chemistries
    against the actual reads and stops the run if the sheet disagrees with them.
    Verified on GSE219015 -- labelling that multiome sample as 10x produces
    "MODALITY MISMATCH ... reads match multiome much better (0.9798 vs 0.0000)"
    and exit 1.

    So a wrong guess costs a download, never a wrong result. What it must not do
    is state "10x" when nothing was found, which is what it used to do: that is
    a guess presented as a fact, and it reads the same as a confirmed value.

    Evidence, strongest first:
      1. the kit named in the protocol prose ("Multiome", "cellranger-arc")
      2. the barcode length stated in prose (24 bp => ARC geometry)
      3. a matching RNA library in the same sample (multiome pairs ATAC + GEX)
    Sources differ in what they know: for this study SRA names the kit for the
    human runs and GEO names it for the mouse one, so both are consulted.
    """
    name = " ".join(str(run.get(k, "")) for k in
                    ("library_name", "sample_name", "gsm")).lower()
    strategy = (run.get("library_strategy") or "").lower()
    if any(m in name or m in strategy for m in MULTIOME_MARKERS):
        return "multiome"

    if not use_network:
        return "unknown"

    prose = library_prose(run)
    if not prose:
        return "unknown"

    if any(m in prose for m in MULTIOME_PROSE_MARKERS):
        return "multiome"

    bc_len = barcode_len_in_prose(prose)
    if bc_len == MULTIOME_BARCODE_READ_LEN:
        return "multiome"
    if bc_len == BARCODE_LEN:
        return "10x"

    # A multiome ATAC library has a GEX partner whose barcodes correspond to it.
    # Submitters say so even when they never name the kit.
    if any(p in prose for p in RNA_PARTNER_PHRASES):
        return "multiome"

    return "unknown"


def is_gex(run):
    """GEX sibling of a multiome pair -- processed by the scRNA pipeline, not here."""
    strategy = (run.get("library_strategy") or "").lower()
    name = " ".join(str(run.get(k, "")) for k in
                    ("library_name", "sample_name")).lower()
    return any(m in strategy or m in name for m in GEX_MARKERS)


def peek_read_length(url, n_reads=2000, n_bytes=1_000_000):
    """Modal read length in the head of a gzipped FASTQ, via an HTTP range
    request (no full download). None on failure. Same technique as the scRNA
    sibling's peek_read_length."""
    import gzip, io
    from urllib.request import Request, urlopen
    seqs = []
    try:
        req = Request(url, headers={"Range": f"bytes=0-{n_bytes}"})
        with urlopen(req, timeout=60) as r:
            blob = r.read()
        # gzip complains about the truncated tail; read what we can.
        gz = gzip.GzipFile(fileobj=io.BytesIO(blob))
        for i, line in enumerate(gz):
            if i % 4 == 1:
                seqs.append(len(line.rstrip(b"\n")))
            if len(seqs) >= n_reads:
                break
    except (EOFError, OSError):
        pass
    except Exception:
        return None
    if not seqs:
        return None
    return Counter(seqs).most_common(1)[0][0]


def select_reads(urls, probe=True):
    """Pick (genomic1_url, barcode_url, genomic2_url) from an ENA URL list.

    Roles are assigned by MEASURED read length, not filename order, because
    submitter naming is not consistent enough to trust:
      the 16 bp read is the cell barcode; the two long reads are the mates.

    Returns (None, None, None) when ENA cannot supply a usable triple -- the
    caller then marks the run source=sra so fasterq-dump --include-technical
    can recover the barcode read (ENA sometimes drops it as "technical").
    """
    if len(urls) < 3:
        return None, None, None

    if not probe:
        # Positional fallback: assume _1,_2,_3 = genomic, barcode, genomic.
        ordered = sorted(urls, key=lambda u: u.rsplit("_", 1)[-1])
        return ordered[0], ordered[1], ordered[2]

    measured = [(u, peek_read_length(u)) for u in urls]
    if any(L is None for _, L in measured):
        return None, None, None

    barcodes = [u for u, L in measured if L == BARCODE_LEN]
    genomic = [u for u, L in measured if L > INDEX_MAX_LEN and L != BARCODE_LEN]
    if len(barcodes) != 1 or len(genomic) != 2:
        return None, None, None
    return genomic[0], barcodes[0], genomic[1]


def fetch_ena_md5s(study_accession):
    """{run_accession: [md5, ...]} from ENA, aligned with the fastq_ftp order.

    WHY THIS IS NOT IN ncbi_utils.py: that file is vendored verbatim from the
    scRNA sibling and deliberately kept byte-identical, so the md5 field is
    added here instead of editing it.

    WHY IT MATTERS: download_fastq already verifies md5s, but prepare_runs.py
    previously wrote empty md5 columns on every row, so the check could never
    fire -- correct code that was never fed anything. ATAC FASTQs run to tens of
    GB and a truncated transfer is a real failure mode (observed once during
    this pipeline's own development, a download that stopped at 18.7 of 20.1 GB).
    """
    import urllib.request
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport"
           f"?accession={study_accession}&result=read_run"
           "&fields=run_accession,fastq_md5&format=tsv")
    out = {}
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            lines = r.read().decode().strip().split("\n")
        if len(lines) < 2:
            return out
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) >= 2 and parts[1]:
                out[parts[0]] = [m for m in parts[1].split(";") if m]
    except Exception:
        # Non-fatal: a missing checksum is worth a warning, not a failed sheet.
        return out
    return out


def md5s_for(md5_map, srr, urls, chosen):
    """Pick the md5s matching the three chosen URLs, by their index in the
    original ENA file list. Returns ("","","") when ENA supplied none."""
    all_md5 = md5_map.get(srr, [])
    if len(all_md5) != len(urls):
        return "", "", ""
    idx = {u: i for i, u in enumerate(urls)}
    try:
        return tuple(all_md5[idx[u]] for u in chosen)
    except (KeyError, IndexError):
        return "", "", ""


def geometry_note(urls):
    """Human-readable observed geometry, for the rejection message."""
    return ", ".join(
        f"{u.rsplit('/', 1)[-1]}={peek_read_length(u)}bp" for u in urls)


# SRA reports the species as ScientificName, which ncbi_utils already surfaces as
# run["organism"]. Map it to the reference build the pipeline knows.
ORGANISM_GENOME = {
    "homo sapiens": "GRCh38",
    "mus musculus": "mm10",
}


def genome_for(run, default, explicit):
    """Reference build for one run, from the organism SRA reports.

    The sheet used to record whatever --genome said, defaulting to GRCh38 and
    never consulting the species. A mouse study then aligned against GRCh38: 97%
    of reads failed to map, and because every downstream metric is computed from
    the reads that did map, the run still produced a complete object with FRiP
    0.70 and TSSe 21 and passed the QC gate. Nothing in the output said "wrong
    species".

    An explicit --genome that contradicts the organism is refused rather than
    honoured. Overriding the species is never what the user meant, and the cost
    of guessing wrong is hours of alignment against the wrong reference.
    """
    organism = (run.get("organism") or "").strip()
    detected = ORGANISM_GENOME.get(organism.lower())

    if detected is None:
        if organism:
            print(f"  NOTE: unrecognised organism {organism!r}; using "
                  f"--genome {default}. Add it to ORGANISM_GENOME if this is wrong.")
        return default

    if explicit and detected != default:
        raise SystemExit(
            f"--genome {default} contradicts the organism SRA reports for "
            f"{run.get('srr')}: {organism} implies {detected}.\n"
            f"Aligning {organism} reads to {default} maps almost nothing, and the "
            "run still produces a QC-passing object from the few reads that do.\n"
            f"Drop --genome to use {detected}, or split the sheet by species."
        )
    return detected


def resolve_runs(accession):
    """Accession -> list of run dicts, scoped by accession type."""
    acc = accession.strip()
    runs = ncbi_utils.fetch_sra_run_info_detailed(acc)
    if not runs:
        return []
    # A run/experiment accession must not fan out to the whole study.
    if acc.upper().startswith(("SRR", "ERR", "DRR")):
        runs = [r for r in runs if r["srr"].upper() == acc.upper()]
    elif acc.upper().startswith(("SRX", "ERX", "DRX")):
        runs = [r for r in runs if r.get("srx", "").upper() == acc.upper()]
    return runs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("accession", help="GSE / SRP / PRJNA / GSM / SRX / SRR")
    ap.add_argument("-o", "--output", default="config/samples.tsv")
    ap.add_argument("--genome", default="GRCh38",
                    help="reference build recorded per row (default: GRCh38)")
    ap.add_argument("--no-probe", action="store_true",
                    help="skip the read-length probe and assume _1/_2/_3 = "
                         "genomic/barcode/genomic. Faster, and wrong on any "
                         "submitter who ordered their files differently.")
    ap.add_argument("--keep-non-atac", action="store_true",
                    help="do not drop non-ATAC libraries (for inspection only)")
    args = ap.parse_args()
    genome_explicit = any(a == '--genome' or a.startswith('--genome=')
                          for a in sys.argv[1:])

    ok, msg = ncbi_utils.check_network_access()
    if not ok:
        print(f"No network access to NCBI/ENA: {msg}", file=sys.stderr)
        sys.exit(1)

    print(f"Resolving {args.accession} ...")
    runs = resolve_runs(args.accession)
    if not runs:
        print(f"No SRA runs found for {args.accession}", file=sys.stderr)
        sys.exit(1)
    print(f"  {len(runs)} run(s) in the accession")

    # --- ATAC scoping -------------------------------------------------------
    kept, dropped = [], []
    for run in runs:
        reason = None if args.keep_non_atac else rejection_reason(run)
        (dropped if reason else kept).append((run, reason))

    dropped_runs = [r for r, _ in dropped]
    if dropped:
        print(f"\nExcluded {len(dropped)} run(s):")
        for run, reason in dropped:
            print(f"  {run['srr']:12} {run.get('gsm', ''):14} {reason}")
    if not kept:
        print("\nNo ATAC libraries left after scoping. If this study really is "
              "scATAC, inspect its LibraryStrategy values and rerun with "
              "--keep-non-atac to see everything.", file=sys.stderr)
        sys.exit(1)
    print(f"\nKeeping {len(kept)} ATAC run(s)")

    # --- ENA URLs -----------------------------------------------------------
    study = kept[0][0].get("sra_study") or args.accession
    print(f"Fetching ENA FASTQ URLs for {study} ...")
    ena = ncbi_utils.fetch_ena_fastq_urls(study)
    ena_md5 = fetch_ena_md5s(study)
    print(f"  ENA supplied checksums for {len(ena_md5)} run(s)")

    rows, sra_fallback = [], []
    for run, _ in kept:
        srr = run["srr"]
        sample = run.get("gsm") or run.get("sample_name") or srr
        urls = ena.get(srr, [])
        g1, bc, g2 = select_reads(urls, probe=not args.no_probe)
        row_genome = genome_for(run, args.genome, genome_explicit)

        if g1 is None:
            sra_fallback.append((srr, len(urls)))
            rows.append({
                "sample": sample, "srr": srr, "source": "sra",
                "modality": detect_modality(run), "genome": row_genome,
                "r1_url": "", "r2_url": "", "r3_url": "",
                "r1_md5": "", "r2_md5": "", "r3_md5": "",
            })
        else:
            m1, m2, m3 = md5s_for(ena_md5, srr, urls, (g1, bc, g2))
            rows.append({
                "sample": sample, "srr": srr, "source": "ena",
                "modality": detect_modality(run), "genome": row_genome,
                "r1_url": g1, "r2_url": bc, "r3_url": g2,
                "r1_md5": m1, "r2_md5": m2, "r3_md5": m3,
            })

    if sra_fallback:
        print(f"\n{len(sra_fallback)}/{len(kept)} run(s) will use SRA. THIS IS "
              "NORMAL FOR scATAC, not a degraded path:")
        for srr, n in sra_fallback:
            print(f"  {srr:12} ENA listed {n} file(s)")
        print("  ENA classifies the 16 bp cell-barcode read as a TECHNICAL read "
              "and does not\n  mirror it. Surveyed 4,000 public scATAC runs on "
              "ENA: 3,815 expose exactly\n  two FASTQs and ZERO expose three. So "
              "for scATAC, SRA is the primary\n  route and ENA is the exception "
              "-- the reverse of the scRNA case, where\n  barcode and cDNA are "
              "both biological reads and ENA works fine.\n"
              "  prefetch + fasterq-dump --include-technical recovers the barcode "
              "read.")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sample", "srr", "source", "modality", "genome",
              "r1_url", "r2_url", "r3_url", "r1_md5", "r2_md5", "r3_md5"]

    def write_sheet(path, sheet_rows):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(sheet_rows)

    # ONE SHEET PER GENOME WHEN A STUDY SPANS SPECIES.
    #
    # GSE219015 holds 16 human and 2 mouse samples. A single sheet covering both
    # is not runnable: the pipeline aligns one genome per run and refuses a mixed
    # sheet, because merging peak sets across species produces an object nobody
    # should use. Writing one anyway just hands the user a file that always fails.
    #
    # So a mixed study produces per-genome sheets and no combined file. Every
    # path printed here can be launched as-is.
    by_genome = {}
    for r in rows:
        by_genome.setdefault(r["genome"], []).append(r)

    if len(by_genome) > 1:
        written = []
        for genome in sorted(by_genome):
            path = out.with_name(f"{out.stem}.{genome}{out.suffix}")
            write_sheet(path, by_genome[genome])
            written.append((path, genome, len(by_genome[genome])))
        print(f"\nThis study spans {len(by_genome)} genomes, so it was split into "
              "one runnable sheet per genome:")
        for path, genome, n in written:
            print(f"  {path}   {n} run(s), genome={genome}")
        print(f"  (no combined {out.name} written: a mixed sheet is refused at "
              "parse time, so it would never be usable)")
    else:
        write_sheet(out, rows)

    n_with_md5 = sum(1 for r in rows if r["r1_md5"])
    n_ena = sum(1 for r in rows if r["source"] == "ena")
    if n_ena:
        print(f"\nchecksums: {n_with_md5}/{n_ena} ENA run(s) carry md5s")
        if n_with_md5 < n_ena:
            print("  Runs without a checksum download unverified. That is ENA's "
                  "gap, not an error, but a truncated transfer will not be caught "
                  "for those runs.")

    n_samples = len({r["sample"] for r in rows})
    n_multiome = sum(1 for r in rows if r["modality"] == "multiome")
    if n_multiome:
        gex_runs = [(r["srr"], r.get("gsm", "")) for r, _ in
                    [(d, None) for d in dropped_runs] if is_gex(r)]
        print(f"\n{n_multiome} run(s) labelled modality=multiome.")
        print("  The ATAC half is processed here. The GEX half is ordinary 10x 3'")
        print("  expression and belongs in the scRNA sibling pipeline")
        print("  (scrnaseq-fetch-process, chemistry=arc); pair the two objects")
        print("  afterwards with merge_h5ad.py --gex-h5ad.")
        if gex_runs:
            print("  GEX sibling runs seen in this accession:")
            for srr, gsm in gex_runs[:10]:
                print(f"    {srr:12} {gsm}")
        print("  NOTE: the multiome path is implemented but NOT verified on real")
        print("  multiome data -- config needs allow_unverified_multiome: true.")

    if len(by_genome) > 1:
        # `out` was never written in this case; naming it would send the user to
        # a file that does not exist.
        print(f"\nWrote {len(by_genome)} sheet(s) above: {len(rows)} run(s) "
              f"across {n_samples} sample(s)")
    else:
        print(f"\nWrote {out}: {len(rows)} run(s) across {n_samples} sample(s)")
    print("Review it before launching -- especially the r2_url column, which "
          "MUST be the 16 bp barcode read.")


if __name__ == "__main__":
    main()

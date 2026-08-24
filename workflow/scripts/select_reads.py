#!/usr/bin/env python3
"""Assign genomic/barcode roles to a run's FASTQ files BY READ LENGTH.

Why not by filename: submitters deposit 10x scATAC under inconsistent names --
R1/R2/R3, R1/R2/I2, R1/I1/R2 -- and SRA's fasterq-dump just numbers them _1.._4.
The one thing that is stable across all of them is the geometry:

    two ~50 bp genomic reads  +  one 16 bp cell-barcode read
    (and sometimes an 8 bp i7 sample index, which we discard)

So we measure the first read of each file and assign roles from that. A file
whose modal read length is 16 is the barcode; the two longest are the genomic
mates, in their original order.

Used by the SRA fallback path in rule download_fastq. The ENA path gets its
per-role URLs from prepare_runs.py, which applies this same rule upstream.
"""
import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

# Accepted cell-barcode read lengths, by chemistry:
#   16 bp = plain 10x single-cell ATAC (barcode occupies the whole read)
#   24 bp = 10x MULTIOME ATAC -- an 8 bp spacer followed by the 16 bp barcode.
#           This is the same fact encoded as MODALITY_OFFSET["multiome"] = 8 in
#           detect_barcode_orientation.py, which slices [8:24] to get the barcode.
#
# THIS LIST EXISTS BECAUSE OF A REAL BUG: it used to be a single BARCODE_LEN = 16,
# so on a multiome run (GSE219015) classify() found zero 16 bp reads, called all
# three reads genomic, and raised -- AFTER a 25-minute download and ~40 minutes of
# extraction, five times over. The two scripts disagreed about what a multiome
# barcode read looks like, and every test fixture used the 10x geometry, so
# nothing caught it.
BARCODE_LENS = (16, 24)
INDEX_MAX_LEN = 10   # i7/i5 sample indexes are 8-10 bp; never a cell barcode


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def modal_read_length(path, n=1000):
    """Modal sequence length over the first n records."""
    lengths = {}
    with _open(path) as fh:
        for i, line in enumerate(fh):
            if i >= n * 4:
                break
            if i % 4 == 1:
                L = len(line.rstrip("\n"))
                lengths[L] = lengths.get(L, 0) + 1
    if not lengths:
        raise ValueError(f"{path} contains no reads")
    return max(lengths, key=lengths.get)


def classify(paths):
    """-> (genomic1, barcode, genomic2). Raises with a readable diagnosis."""
    measured = [(p, modal_read_length(p)) for p in paths]
    shown = ", ".join(f"{Path(p).name}={L}bp" for p, L in measured)

    barcodes = [p for p, L in measured if L in BARCODE_LENS]
    genomic = [p for p, L in measured
               if L > INDEX_MAX_LEN and L not in BARCODE_LENS]

    if len(barcodes) != 1:
        raise ValueError(
            f"expected exactly one cell-barcode read of {BARCODE_LENS[0]}bp "
            f"(10x) or {BARCODE_LENS[1]}bp (multiome), found {len(barcodes)}. "
            f"Observed geometry: {shown}.\n"
            "This does not look like 10x single-cell ATAC or Multiome ATAC. If "
            "it is a different chemistry, this pipeline does not support it."
        )
    if len(genomic) != 2:
        raise ValueError(
            f"expected exactly two genomic reads, found {len(genomic)}. "
            f"Observed geometry: {shown}."
        )
    return genomic[0], barcodes[0], genomic[1]


def place(src, dst, compress):
    """Put the chosen read where the rule expects it.

    compress=False (the default, used when keep_fastq is false):
        MOVE the file. No compression at all. The FASTQs are deleted by
        Snakemake temp() as soon as chromap has read them, and chromap reads
        plain FASTQ fine -- so compressing them costs hours to produce something
        thrown away minutes later. Measured: Python's gzip does ~25 MB/s, so
        ~460 GB took roughly 5 hours, longer than the download and the alignment
        put together.

        A move is also near-instant when src and dst share a filesystem, and
        falls back to a copy when they do not.

    compress=True (keep_fastq: true):
        gzip it, because the file is meant to persist. Uses pigz when available
        (parallel, ~200 MB/s on these nodes) and falls back to Python gzip.
    """
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)

    if not compress:
        if str(src).endswith(".gz"):
            # already compressed upstream; nothing to undo
            shutil.move(str(src), str(dst))
        else:
            shutil.move(str(src), str(dst))
        return

    if str(src).endswith(".gz"):
        shutil.copyfile(src, dst)
        return
    pigz = shutil.which("pigz")
    if pigz:
        with open(dst, "wb") as fo:
            subprocess.run([pigz, "-p", "8", "-c", str(src)], stdout=fo, check=True)
    else:
        with open(src, "rb") as fi, gzip.open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--from-dir", required=True,
                    help="directory of FASTQs produced by fasterq-dump")
    ap.add_argument("--r1", required=True, help="output: genomic read 1 (.gz)")
    ap.add_argument("--r2", required=True, help="output: cell barcode read (.gz)")
    ap.add_argument("--r3", required=True, help="output: genomic read 2 (.gz)")
    ap.add_argument("--compress", action="store_true",
                    help="gzip the outputs (pigz when available). Off by "
                         "default: with keep_fastq=false the FASTQs are deleted "
                         "right after chromap reads them, so compressing them "
                         "is hours of work thrown away.")
    args = ap.parse_args()

    found = sorted(Path(args.from_dir).glob("*.fastq")) + \
            sorted(Path(args.from_dir).glob("*.fastq.gz"))
    if not found:
        print(f"No FASTQ files in {args.from_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        g1, bc, g2 = classify(found)
    except ValueError as e:
        print(f"Read-geometry check failed in {args.from_dir}: {e}", file=sys.stderr)
        sys.exit(1)

    bc_len = modal_read_length(bc)
    chem = "10x (16bp barcode)" if bc_len == 16 else "multiome (24bp: 8bp spacer + 16bp barcode)"
    print(f"genomic1={Path(g1).name}  barcode={Path(bc).name}  genomic2={Path(g2).name}")
    print(f"barcode read is {bc_len}bp -> {chem}")
    print("  (cross-check this against the `modality` column in samples.tsv; "
          "detect_barcode_orientation will also probe both chemistries)")
    place(g1, args.r1, args.compress)
    place(bc, args.r2, args.compress)
    place(g2, args.r3, args.compress)
    print("compressed with pigz/gzip" if args.compress
          else "moved uncompressed (deleted by temp() after chromap)")


if __name__ == "__main__":
    main()

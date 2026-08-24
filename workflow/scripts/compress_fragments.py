#!/usr/bin/env python3
"""Sort, bgzip and tabix-index a chromap fragment BED.

chromap writes an unsorted plain-text BED. snapATAC2, IGV and any random-access
consumer need it coordinate-sorted, BGZF-compressed and tabix-indexed.

We do this through pysam rather than the bgzip/tabix CLIs, which keeps htslib a
Python dependency instead of a second binary to install. scATACpipe emits plain
gzip here with the note that "bgzip creates another layer of complexity for
Python to handle" -- that leaves the fragment file unindexable, so we do not
copy it.

Sorting is done with the system `sort` (external merge sort, bounded memory)
rather than in Python, because a deep ATAC library can exceed a comfortable
in-memory sort.
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw", required=True, help="unsorted fragment BED from chromap")
    ap.add_argument("--out", required=True, help="output fragments.tsv.gz (bgzipped)")
    ap.add_argument("--sort-mem", default="4G", help="memory budget for `sort -S`")
    args = ap.parse_args()

    raw, out = Path(args.raw), Path(args.out)
    if not raw.exists() or raw.stat().st_size == 0:
        print(f"chromap produced no fragments: {raw} is missing or empty.\n"
              "Common causes: the barcode whitelist orientation is wrong (check "
              "qc/barcodes/<sample>.orientation.txt), or the genomic reads were "
              "passed in the barcode slot.", file=sys.stderr)
        sys.exit(1)

    out.parent.mkdir(parents=True, exist_ok=True)
    sorted_bed = out.with_suffix("").with_suffix(".sorted.bed")

    # LC_ALL=C makes the sort byte-wise and locale-independent, which is what
    # tabix expects.
    env = dict(os.environ, LC_ALL="C")
    with open(sorted_bed, "w") as fh:
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", "-S", args.sort_mem, str(raw)],
            stdout=fh, env=env, check=True,
        )

    pysam.tabix_compress(str(sorted_bed), str(out), force=True)
    pysam.tabix_index(str(out), preset="bed", force=True)
    sorted_bed.unlink()

    # Prove the index actually works rather than trusting that it was written.
    with pysam.TabixFile(str(out)) as tbx:
        contigs = list(tbx.contigs)
    if not contigs:
        print(f"tabix index for {out} contains no contigs -- the fragment file "
              "is malformed.", file=sys.stderr)
        sys.exit(1)

    n_bytes = out.stat().st_size
    print(f"wrote {out} ({n_bytes:,} bytes, {len(contigs)} contigs) + .tbi")


if __name__ == "__main__":
    main()

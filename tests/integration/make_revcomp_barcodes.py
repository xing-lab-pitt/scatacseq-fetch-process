#!/usr/bin/env python3
"""Reverse-complement a barcode FASTQ, to exercise the revcomp code path.

WHY: whether a 10x barcode whitelist matches the reads as written or
reverse-complemented depends on the SEQUENCER, not the assay. Both datasets run
through this pipeline so far chose 'forward' (match rates 0.939 and 0.899, with
revcomp at 0.000), so the revcomp branch has never actually executed on real
data -- even though it certainly occurs in the wild (ENCODE_scatac ships the
reverse-complemented whitelist as its 10x default).

This flips a known-good dataset in silico so the branch runs. The pipeline should
then (a) detect revcomp, and (b) produce the SAME biological result as the
unflipped run. (b) is the real assertion: detecting the flip is useless if the
resulting whitelist does not actually let chromap correct barcodes.

Quality strings are reversed alongside the sequence, so per-base qualities stay
attached to their bases -- otherwise FastQC metrics would shift for a reason
unrelated to what is being tested.
"""
import argparse
import gzip
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-fastq", required=True)
    ap.add_argument("--out-fastq", required=True)
    args = ap.parse_args()

    n = 0
    Path(args.out_fastq).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.in_fastq, "rt") as fi, gzip.open(args.out_fastq, "wt") as fo:
        while True:
            head = fi.readline()
            if not head:
                break
            seq = fi.readline().rstrip("\n")
            plus = fi.readline()
            qual = fi.readline().rstrip("\n")
            fo.write(head)
            fo.write(revcomp(seq) + "\n")
            fo.write(plus)
            fo.write(qual[::-1] + "\n")   # reverse quality to match the sequence
            n += 1

    print(f"reverse-complemented {n:,} barcode reads -> {args.out_fastq}")


if __name__ == "__main__":
    main()

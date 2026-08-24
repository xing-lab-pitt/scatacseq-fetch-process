#!/usr/bin/env python3
"""Filter a MACS3 peak set: drop blacklisted regions and non-primary contigs.

WHY THIS IS A SEPARATE RULE, not folded into macs3_peaks:
MACS3's raw output is kept on disk so the drop can be audited -- "35,037 called,
34,「N」 kept" is a QC signal in its own right, and an unexpectedly large drop
means something is wrong with the peak call, not with the filter. Keeping them
separate also lets Snakemake re-run just the filter when the blacklist changes,
without re-calling peaks.

TWO FILTERS, both optional and independently controlled:

1. BLACKLIST (the one that matters). The ENCODE blacklist (Amemiya et al., Sci
   Rep 2019) marks regions that produce high signal in *any* ChIP/ATAC
   experiment regardless of biology -- collapsed repeats, centromeric satellite,
   and low-mappability regions. Peaks there are artifacts: they attract reads
   from everywhere in the genome, so they look like the strongest peaks in the
   sample and pull FRiP upward while carrying no signal. hg38 v2 covers ~227 Mb
   across 636 regions.

2. PRIMARY CONTIGS. Unplaced scaffolds (GL000194.1 etc.) and chrM. Small in
   number but they are not interpretable as regulatory elements, and chrM in
   particular is an accessibility artifact of mitochondrial copy number rather
   than chromatin state.

CONTIG NAMING is checked explicitly. A blacklist using `chr1` against peaks
using `1` would silently filter nothing and quietly leave the artifacts in --
the exact failure mode this script exists to prevent, so it refuses instead.
"""
import argparse
import re
import sys
from pathlib import Path

import pandas as pd

NARROWPEAK_COLS = ["chrom", "start", "end", "name", "score",
                   "strand", "signal", "pvalue", "qvalue", "summit"]
# chr1..chr22, chrX, chrY -- deliberately excludes chrM and unplaced scaffolds.
PRIMARY_RE = re.compile(r"^chr([0-9]{1,2}|X|Y)$")


def read_narrowpeak(path):
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    df.columns = NARROWPEAK_COLS[:df.shape[1]]
    return df


def read_bed(path):
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"], comment="#")
    return df


def overlaps_blacklist(peaks, black):
    """Boolean mask: True where a peak overlaps any blacklist region.

    Plain interval sweep per chromosome rather than pyranges: the blacklist is
    636 rows and peaks are ~10^5, so this is instant, and it keeps the script
    readable and dependency-light.
    """
    mask = pd.Series(False, index=peaks.index)
    for chrom, bl in black.groupby("chrom"):
        idx = peaks.index[peaks["chrom"] == chrom]
        if len(idx) == 0:
            continue
        starts = bl["start"].to_numpy()
        ends = bl["end"].to_numpy()
        p_start = peaks.loc[idx, "start"].to_numpy()
        p_end = peaks.loc[idx, "end"].to_numpy()
        # A peak overlaps if it starts before some blacklist region ends AND
        # ends after that region starts.
        hit = ((p_start[:, None] < ends[None, :]) &
               (p_end[:, None] > starts[None, :])).any(axis=1)
        mask.loc[idx] = hit
    return mask


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--peaks", required=True, help="MACS3 narrowPeak input")
    ap.add_argument("--out", required=True, help="filtered narrowPeak output")
    ap.add_argument("--blacklist", default="",
                    help="BED of blacklisted regions; empty disables the filter")
    ap.add_argument("--genome-fai", default="",
                    help="samtools .fai for the genome the peaks were called on. "
                         "Used to verify the blacklist belongs to this genome.")
    ap.add_argument("--primary-chroms-only", action="store_true",
                    help="keep only chr1-22, chrX, chrY (drops chrM + scaffolds)")
    ap.add_argument("--report", default="", help="optional TSV summary")
    args = ap.parse_args()

    peaks = read_narrowpeak(args.peaks)
    n_start = len(peaks)
    dropped = {}

    if args.primary_chroms_only:
        keep = peaks["chrom"].astype(str).str.match(PRIMARY_RE)
        dropped["non_primary_contig"] = int((~keep).sum())
        peaks = peaks[keep]

    if args.blacklist:
        bl_path = Path(args.blacklist)
        if not bl_path.exists():
            print(f"Blacklist not found: {bl_path}", file=sys.stderr)
            sys.exit(1)
        black = read_bed(bl_path)

        # Refuse on a naming mismatch rather than filtering nothing silently.
        shared = set(peaks["chrom"].astype(str)) & set(black["chrom"].astype(str))
        if not shared:
            print(
                f"Contig naming mismatch: no chromosome name is shared between\n"
                f"  peaks     ({sorted(set(peaks['chrom'].astype(str)))[:4]} ...)\n"
                f"  blacklist ({sorted(set(black['chrom'].astype(str)))[:4]} ...)\n"
                "A blacklist that matches nothing would silently leave every "
                "artifact peak in place, so this is treated as a hard error. Use "
                "a blacklist built for the same reference naming convention "
                "(hg38-blacklist.v2 uses the 'chr' prefix).",
                file=sys.stderr)
            sys.exit(1)

        # Sharing SOME contig names is not enough. hg38 and mm10 both use
        # chr1-chr19, chrX and chrY, so a human blacklist against mouse peaks
        # passes the check above, then drops the wrong ~9,000 peaks while leaving
        # the real artifacts in place. Compare against the genome's own contig
        # list instead: a blacklist naming chr20 is not an mm10 blacklist.
        if args.genome_fai and Path(args.genome_fai).exists():
            genome_chroms = {l.split("\t")[0]
                             for l in Path(args.genome_fai).read_text().splitlines()
                             if l.strip()}
            alien = sorted(set(black["chrom"].astype(str)) - genome_chroms)
            if alien:
                print(
                    f"Blacklist does not belong to this genome.\n"
                    f"  blacklist : {bl_path}\n"
                    f"  genome    : {args.genome_fai}\n"
                    f"  {len(alien)} blacklist contig(s) absent from the genome, "
                    f"e.g. {alien[:5]}\n"
                    "Sharing chr1-chr19/chrX/chrY is not evidence of a match -- "
                    "human and mouse share those names. Point "
                    "references.<genome>.blacklist at the right file.",
                    file=sys.stderr)
                sys.exit(1)

        hit = overlaps_blacklist(peaks, black)
        dropped["blacklisted"] = int(hit.sum())
        peaks = peaks[~hit]

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    peaks.to_csv(args.out, sep="\t", header=False, index=False)

    n_end = len(peaks)
    pct = 100.0 * (n_start - n_end) / n_start if n_start else 0.0
    print(f"peaks: {n_start:,} called -> {n_end:,} kept "
          f"({n_start - n_end:,} dropped, {pct:.2f}%)")
    for reason, n in dropped.items():
        print(f"  {reason}: {n:,}")

    if args.report:
        rows = [{"metric": "peaks_called", "value": n_start},
                {"metric": "peaks_kept", "value": n_end}]
        rows += [{"metric": f"dropped_{k}", "value": v} for k, v in dropped.items()]
        Path(args.report).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(args.report, sep="\t", index=False)

    if n_end == 0:
        print("\nAll peaks were filtered out. Check the blacklist is the right "
              "genome build.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

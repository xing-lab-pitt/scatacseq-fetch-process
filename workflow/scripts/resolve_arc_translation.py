#!/usr/bin/env python3
"""Build the 10x Multiome (ARC) ATAC->GEX barcode translation table.

WHY A TABLE IS NECESSARY (verified, not assumed):
In 10x Multiome one gel bead carries two different barcode sequences, one used by
the ATAC library and one by the GEX library. Measured on the real ARC whitelists:

    ATAC list and GEX list      736,319 entries each  (equal -- pairing is by line)
    shared barcodes             0 of 736,319          (disjoint universes)
    is ATAC == revcomp(GEX)?    no                    (unrelated sequences)

So the mapping is arbitrary: it exists ONLY as line correspondence between the
two files and cannot be computed from a barcode's sequence. Hence a lookup table,
handed to chromap via --barcode-translate, which makes chromap emit fragments
keyed by the GEX barcode so ATAC and GEX objects share obs_names.

THE TRAP THIS GUARDS AGAINST:
Cell Ranger ARC ships two DIFFERENT files with the same basename
737K-arc-v1.txt -- lib/python/atac/barcodes/ (ATAC) and
lib/python/cellranger/barcodes/ (GEX). Silently passing the GEX list where the
ATAC one belongs produces a run that completes and looks fine but is keyed
wrongly throughout. The checks below are here to make that impossible.

Built once and cached beside the whitelists, same pattern as the chromap index.
"""
import argparse
import gzip
import sys
from pathlib import Path

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_whitelist(path):
    with _open(path) as fh:
        return [line.strip().split()[0] for line in fh if line.strip()]


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def validate(atac, gex, atac_path, gex_path):
    """Every way these two files can be wrong, checked explicitly."""
    problems = []

    if len(atac) != len(gex):
        problems.append(
            f"line counts differ: {len(atac):,} (ATAC) vs {len(gex):,} (GEX). "
            "The pairing is BY LINE NUMBER, so unequal files mean these are not "
            "a matched ARC pair and every translation would be wrong.")
        return problems  # nothing else is meaningful if lengths differ

    if atac == gex:
        problems.append(
            "the two whitelists are IDENTICAL. You have passed the same file "
            "twice. Cell Ranger ARC ships two different files both named "
            "737K-arc-v1.txt -- take the ATAC one from lib/python/atac/barcodes/ "
            "and the GEX one from lib/python/cellranger/barcodes/.")

    shared = set(atac) & set(gex)
    frac = len(shared) / len(atac) if atac else 0
    if frac > 0.01:
        problems.append(
            f"{len(shared):,} barcodes ({frac:.1%}) appear in BOTH lists. Real "
            "ARC whitelists are disjoint universes (measured: 0 shared). A large "
            "overlap means at least one file is not what it claims to be.")

    bad_len = {len(b) for b in atac[:1000]} | {len(b) for b in gex[:1000]}
    if bad_len != {16}:
        problems.append(f"expected 16 bp barcodes, saw lengths {sorted(bad_len)}")

    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--atac-whitelist", required=True,
                    help="ARC *ATAC* whitelist (lib/python/atac/barcodes/)")
    ap.add_argument("--gex-whitelist", required=True,
                    help="ARC *GEX* whitelist (lib/python/cellranger/barcodes/)")
    ap.add_argument("-o", "--output", required=True,
                    help="two-column TSV: atac_barcode <TAB> gex_barcode")
    args = ap.parse_args()

    for p in (args.atac_whitelist, args.gex_whitelist):
        if not Path(p).exists():
            print(f"Whitelist not found: {p}", file=sys.stderr)
            sys.exit(1)

    atac = read_whitelist(args.atac_whitelist)
    gex = read_whitelist(args.gex_whitelist)
    print(f"ATAC ARC whitelist: {len(atac):,} barcodes  ({args.atac_whitelist})")
    print(f"GEX  ARC whitelist: {len(gex):,} barcodes  ({args.gex_whitelist})")

    problems = validate(atac, gex, args.atac_whitelist, args.gex_whitelist)
    if problems:
        print("\nARC whitelist pair failed validation:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        sys.exit(1)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    # COLUMN ORDER IS target<TAB>source -- i.e. GEX FIRST, ATAC SECOND.
    #
    # This is counterintuitive and it is NOT documented in `chromap -h`, which
    # says only "Convert barcode to the specified sequences during output".
    # From chromap's source (src/barcode_translator.h, ProcessTranslateFileLine):
    #
    #     to   = line.substr(0, i);          // first column  = TARGET
    #     from = line.substr(i + 1, ...);    // second column = SOURCE (the key)
    #
    # so the second column is what chromap looks up and the first is what it
    # writes out. Getting this backwards does NOT silently mistranslate -- chromap
    # aborts with "Barcode does not exist in the translation table." -- but it
    # cost a full debugging cycle to establish, so it is written down here.
    # Separator may be a tab or a comma.
    with open(out, "w") as fh:
        for a, g in zip(atac, gex):
            fh.write(f"{g}\t{a}\n")

    print(f"\nWrote {out}: {len(atac):,} pairs, written as GEX<TAB>ATAC")
    print(f"  (chromap reads column 2 as the key and emits column 1)")
    print(f"  example: ATAC {atac[0]} will be emitted as GEX {gex[0]}")


if __name__ == "__main__":
    main()

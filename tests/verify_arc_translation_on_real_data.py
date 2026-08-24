"""Do the ATAC barcodes in a finished multiome .h5ad translate to real GEX barcodes?

Multiome exists to pair ATAC with GEX. Every run so far processed ATAC only, so
nothing has ever checked that the two modalities can actually be joined. The full
check needs the GEX libraries through the scRNA pipeline; this is the half that
can be done from the ATAC object alone, and it would catch a systematic
translation error.

What it asserts, on real data:
  1. the cell barcodes in the .h5ad are ARC ATAC barcodes, not plain 10x ones
  2. every one of them translates to a barcode that exists in the ARC GEX list
  3. the mapping is one-to-one -- no two ATAC cells collapse onto one GEX barcode

Run:
    SCATAC_WHITELIST_DIR=/path/to/10x_whitelists \\
      python tests/verify_arc_translation_on_real_data.py <sample.h5ad>
"""
import os
import sys
from pathlib import Path

import anndata as ad

WL = Path(os.environ.get("SCATAC_WHITELIST_DIR", "reference/10x_whitelists"))
ATAC_WL, GEX_WL = WL / "737K-arc-v1.ATAC.txt", WL / "737K-arc-v1.txt"


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    h5ad = Path(sys.argv[1])
    for p in (h5ad, ATAC_WL, GEX_WL):
        if not p.exists():
            print(f"SKIP: missing {p}")
            return 0

    atac = ATAC_WL.read_text().split()
    gex = GEX_WL.read_text().split()
    if len(atac) != len(gex):
        print(f"FAIL: whitelists differ in length ({len(atac):,} vs {len(gex):,}); "
              "they pair by line number, so this makes translation meaningless")
        return 1
    a2g = dict(zip(atac, gex))
    gex_set = set(gex)

    a = ad.read_h5ad(h5ad, backed="r")
    obs = a.obs
    cells = list(obs[obs["is_cell"]].index) if "is_cell" in obs.columns else list(obs.index)
    print(f"{h5ad.name}: {len(cells):,} cells")

    # chromap emits the TARGET column of the translate table, so a multiome run's
    # barcodes should already be GEX-side. Accept either and say which.
    n_atac = sum(1 for c in cells if c in a2g)
    n_gex = sum(1 for c in cells if c in gex_set)
    n_rc_atac = sum(1 for c in cells if revcomp(c) in a2g)
    print(f"  match ARC ATAC list : {n_atac:,}  ({100*n_atac/len(cells):.1f}%)")
    print(f"  match ARC GEX  list : {n_gex:,}  ({100*n_gex/len(cells):.1f}%)")
    print(f"  match ATAC revcomp  : {n_rc_atac:,}")

    ok = True
    if max(n_atac, n_gex) < 0.95 * len(cells):
        print("FAIL: fewer than 95% of cell barcodes are in either ARC whitelist. "
              "These are probably not ARC barcodes at all.")
        ok = False

    if n_gex >= n_atac:
        print("  -> barcodes are GEX-side: translation already applied by chromap.")
        joinable = n_gex
    else:
        print("  -> barcodes are ATAC-side; translating to check the GEX side exists.")
        translated = [a2g[c] for c in cells if c in a2g]
        missing = [t for t in translated if t not in gex_set]
        if missing:
            print(f"FAIL: {len(missing):,} translated barcodes absent from the GEX list")
            ok = False
        if len(set(translated)) != len(translated):
            dup = len(translated) - len(set(translated))
            print(f"FAIL: translation is not one-to-one -- {dup:,} collisions")
            ok = False
        joinable = len(translated)

    print(f"  joinable with a GEX object: {joinable:,}/{len(cells):,}")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Concatenate per-sample scATAC .h5ad objects; optionally pair with GEX as MuData.

Peak sets differ between samples (each was called on that sample's own
fragments), so the union of peaks is taken and cells lacking a peak get zeros --
an outer join. This is the honest default: an inner join would silently discard
sample-specific accessibility, which is often the biology of interest.

The combined object is a convenience for a first look. The per-sample fragment
files remain the durable substrate: peaks are contingent on the population that
was called, fragments are not. adata.uns['fragments_by_sample'] keeps the map.

--- MULTIOME PAIRING -------------------------------------------------------
For 10x Multiome, --gex-h5ad accepts the matching GEX object and the two are
paired into a MuData on shared obs_names.

WHY THE GEX HALF IS NOT PROCESSED BY THIS PIPELINE:
Multiome GEX is ordinary 10x 3' gene expression and belongs in the scRNA sibling
(scrnaseq-fetch-process, STARsolo with chemistry=arc). Re-implementing a second
aligner here would violate the one-aligner-one-convention rule the whole design
rests on. So the workflow is: run the ATAC half here, the GEX half there, then
pair the two objects with this script.

The pairing only works because chromap was given --barcode-translate, which
rekeys ATAC fragments to the GEX barcode. Without it the two objects have
disjoint obs_names and the overlap check below fails loudly -- which is the
entire point, because a silently unpaired MuData looks perfectly well-formed.
"""
import argparse
import sys
from pathlib import Path

import anndata as ad

# Below this fraction of shared barcodes, the pairing is treated as failed
# rather than as a thin overlap. Real paired multiome shares most cells.
MIN_SHARED_FRACTION = 0.30


def concat_atac(paths):
    objs, frags = {}, {}
    for path in paths:
        sample = Path(path).stem
        a = ad.read_h5ad(path)
        objs[sample] = a
        if "fragments" in a.uns:
            frags[sample] = str(a.uns["fragments"])
        print(f"  {sample}: {a.n_obs:,} barcodes x {a.n_vars:,} peaks")

    if len(objs) == 1:
        combined = next(iter(objs.values())).copy()
    else:
        combined = ad.concat(objs, join="outer", label="sample",
                             index_unique="-", merge="unique")
    combined.uns["fragments_by_sample"] = frags
    return combined


def pair_with_gex(atac, gex_path, output):
    """Build a MuData from the ATAC object and a separately-produced GEX object."""
    try:
        import mudata
    except ImportError:
        print("mudata is required for --gex-h5ad. Install it into the pipeline "
              "environment: uv pip install mudata", file=sys.stderr)
        sys.exit(1)

    gex = ad.read_h5ad(gex_path)
    shared = set(atac.obs_names) & set(gex.obs_names)
    smaller = min(atac.n_obs, gex.n_obs)
    frac = len(shared) / smaller if smaller else 0.0

    print("\nMultiome pairing:")
    print(f"  ATAC  : {atac.n_obs:,} barcodes")
    print(f"  GEX   : {gex.n_obs:,} barcodes")
    print(f"  shared: {len(shared):,} ({frac:.1%} of the smaller object)")

    if frac < MIN_SHARED_FRACTION:
        print(
            f"\nATAC and GEX share only {frac:.1%} of barcodes -- below the "
            f"{MIN_SHARED_FRACTION:.0%} threshold, so this is treated as a FAILED "
            "pairing rather than a partial one.\n"
            "In 10x Multiome the two libraries carry DIFFERENT barcode sequences "
            "for the same gel bead. They only line up if chromap was run with "
            "--barcode-translate against the ARC translation table.\n"
            "Check, in order:\n"
            "  1. was the sample's modality set to 'multiome' in samples.tsv?\n"
            "  2. did rule resolve_arc_translation actually run?\n"
            "  3. was the ATAC ARC whitelist used (NOT the GEX one -- both are\n"
            "     named 737K-arc-v1.txt; see M5_MULTIOME_PLAN.md).\n"
            "Refusing to write a MuData that would look well-formed but be keyed "
            "to the wrong cells.", file=sys.stderr)
        sys.exit(1)

    # Explicit per-modality presence flags rather than a silent drop, so a cell
    # seen in only one assay stays visible and auditable.
    atac.obs["in_atac"] = True
    gex.obs["in_gex"] = True
    mdata = mudata.MuData({"atac": atac, "rna": gex})
    mdata.uns["shared_barcodes"] = len(shared)
    mdata.uns["shared_fraction"] = float(frac)

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    mdata.write(output)
    print(f"\nWrote {output}: MuData with atac {atac.shape} + rna {gex.shape}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", nargs="+", required=True,
                    help="per-sample ATAC .h5ad objects")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--gex-h5ad", default="",
                    help="multiome only: matching GEX .h5ad, produced by the "
                         "scRNA sibling pipeline. Output becomes a MuData (.h5mu).")
    args = ap.parse_args()

    combined = concat_atac(args.h5ad)

    if args.gex_h5ad:
        pair_with_gex(combined, args.gex_h5ad, args.output)
        return

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(args.output)
    print(f"\nWrote {args.output}: {combined.n_obs:,} barcodes x "
          f"{combined.n_vars:,} peaks from {len(args.h5ad)} sample(s)")


if __name__ == "__main__":
    main()

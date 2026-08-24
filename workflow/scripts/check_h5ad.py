#!/usr/bin/env python3
"""Structural check on a per-sample scATAC .h5ad.

The QC gate answers "is this sample good data?". This answers the prior question
"is this object even well-formed?" -- a malformed object with plausible numbers
is the failure mode that wastes the most downstream time, because it surfaces as
a confusing error three analyses later.

Every failure names its CAUSE, not just the symptom. Analogue of check_layers.py
in the scRNA sibling, which verifies the Velocyto layers are actually populated.
"""
import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np

REQUIRED_VAR = ["chrom", "start", "end"]
REQUIRED_OBS = ["n_fragments", "frip", "tss_enrichment",
                "nucleosome_signal", "frac_mito", "is_cell"]


def check(path):
    problems = []
    a = ad.read_h5ad(path)

    if a.n_obs == 0:
        problems.append("no barcodes (n_obs=0) -- chromap produced no usable cells")
    if a.n_vars == 0:
        problems.append("no peaks (n_vars=0) -- MACS3 called nothing; check the "
                        "fragment file is non-empty and the genome size is right")

    for col in REQUIRED_VAR:
        if col not in a.var.columns:
            problems.append(f"var is missing '{col}' -- peak coordinates incomplete")

    for col in REQUIRED_OBS:
        if col not in a.obs.columns:
            problems.append(f"obs is missing QC column '{col}'")
            continue
        vals = a.obs[col]
        if vals.isna().all():
            problems.append(f"obs['{col}'] is entirely NaN -- the metric was never "
                            "computed, only allocated")

    if a.X is not None and a.n_obs and a.n_vars:
        nnz = a.X.nnz if hasattr(a.X, "nnz") else int(np.count_nonzero(a.X))
        if nnz == 0:
            problems.append("X has zero non-zero entries -- the peak x cell matrix "
                            "is empty. Usually a chrom-name mismatch between the "
                            "fragments and the peaks (chr1 vs 1).")

    frag = a.uns.get("fragments")
    if not frag:
        problems.append("uns['fragments'] is unset -- the durable fragment file "
                        "path was not recorded")
    elif not Path(frag).exists():
        problems.append(f"uns['fragments'] points at {frag}, which does not exist")

    return a, problems


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("h5ad")
    ap.add_argument("--out", required=True, help="marker file written on success")
    args = ap.parse_args()

    a, problems = check(args.h5ad)

    print(f"{args.h5ad}: {a.n_obs:,} barcodes x {a.n_vars:,} peaks")
    if "is_cell" in a.obs.columns:
        print(f"  called cells: {int(a.obs['is_cell'].sum()):,}")
    for key in ("frac_reads_in_peaks", "median_tss_enrichment",
                "valid_barcode_frac", "duplicate_rate"):
        if key in a.uns:
            print(f"  {key}: {a.uns[key]}")

    if problems:
        print(f"\n{args.h5ad} FAILED {len(problems)} structural check(s):",
              file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        sys.exit(1)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text("ok\n")
    print("structure OK")


if __name__ == "__main__":
    main()

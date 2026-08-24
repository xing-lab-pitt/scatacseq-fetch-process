#!/usr/bin/env python3
"""Deterministic QC gate over per-sample scATAC .h5ad objects.

Applies FIXED thresholds (passed in from config.yaml) to each sample and writes a
pass/fail table plus a list of passing samples. It does NOT choose thresholds
itself, so the decision is reproducible and reviewable in version control --
scATAC-pro's author makes exactly this argument: the cutoffs cannot be
standardised across tissues, so they must at least be explicit.

Metric provenance (each from the cheapest tool that can produce it):
  valid_barcode_frac, duplicate_rate  <- chromap --summary
  tss_enrichment, nucleosome_signal   <- snapATAC2 (peak-independent)
  frac_reads_in_peaks (bulk FRiP)     <- pyranges overlap, needs the peak set
  n_cells                             <- the two-threshold cell call

Always exits 0: it produces the report. Enforcement (stop vs. warn) is a
separate rule so the report is never deleted by a failing gate.
"""
import argparse
import csv
import sys
from pathlib import Path

import anndata as ad

# uns key -> (threshold name, direction). "min" = fail when below, "max" = above.
CHECKS = [
    ("frac_reads_in_peaks",   "min_frac_reads_in_peaks",  "min"),
    ("median_tss_enrichment", "min_tss_enrichment",       "min"),
    ("n_cells",               "min_estimated_cells",      "min"),
    ("valid_barcode_frac",    "min_valid_barcode_frac",   "min"),
    ("duplicate_rate",        "max_duplicate_rate",       "max"),
    # Every other metric here is computed from the fragment file, which contains
    # only reads that mapped. They therefore describe the surviving reads and say
    # nothing about how many were lost. A run of GSE219015 mapped 0.7% of its
    # reads and passed every check above with FRiP 0.70 and TSSe 21.
    ("mapping_rate",          "min_mapping_rate",         "min"),
]


def evaluate(h5ad_path, thresholds):
    a = ad.read_h5ad(h5ad_path, backed="r")
    metrics = {key: a.uns.get(key) for key, _, _ in CHECKS}
    reasons = []

    for key, thr_name, direction in CHECKS:
        thr = thresholds.get(thr_name)
        if thr is None:
            continue
        val = metrics.get(key)
        if val is None:
            reasons.append(f"{key}:missing")
        elif direction == "min" and float(val) < thr:
            reasons.append(f"{key}={float(val):.4g}<{thr}")
        elif direction == "max" and float(val) > thr:
            reasons.append(f"{key}={float(val):.4g}>{thr}")

    return metrics, ("PASS" if not reasons else "FAIL"), ";".join(reasons)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("entries", nargs="+", help="sample=<path/to/sample.h5ad>")
    ap.add_argument("--report", required=True)
    ap.add_argument("--passlist", required=True)
    ap.add_argument("--min-frac-reads-in-peaks", type=float)
    ap.add_argument("--min-tss-enrichment", type=float)
    ap.add_argument("--min-estimated-cells", type=float)
    ap.add_argument("--min-valid-barcode-frac", type=float)
    ap.add_argument("--max-duplicate-rate", type=float)
    ap.add_argument("--min-mapping-rate", type=float)
    args = ap.parse_args()

    thresholds = {
        "min_frac_reads_in_peaks": args.min_frac_reads_in_peaks,
        "min_tss_enrichment": args.min_tss_enrichment,
        "min_estimated_cells": args.min_estimated_cells,
        "min_valid_barcode_frac": args.min_valid_barcode_frac,
        "max_duplicate_rate": args.max_duplicate_rate,
        "min_mapping_rate": args.min_mapping_rate,
    }

    rows, passing = [], []
    for entry in args.entries:
        sample, _, path = entry.partition("=")
        metrics, status, reason = evaluate(path, thresholds)
        row = {"sample": sample}
        row.update({k: metrics.get(k) for k, _, _ in CHECKS})
        row["status"] = status
        row["reasons"] = reason
        rows.append(row)
        if status == "PASS":
            passing.append(sample)

    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    Path(args.passlist).write_text("\n".join(passing) + ("\n" if passing else ""))

    print(f"QC gate: {len(passing)}/{len(rows)} passed. Report: {args.report}")
    for r in rows:
        if r["status"] == "FAIL":
            print(f"  FAIL {r['sample']}: {r['reasons']}", file=sys.stderr)
    sys.exit(0)  # report always succeeds; enforcement is a separate rule


if __name__ == "__main__":
    main()

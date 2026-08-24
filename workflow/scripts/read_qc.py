#!/usr/bin/env python3
"""Flag raw-read quality per run from FastQC output. REPORTING ONLY, never fatal.

Role-aware, because the three reads of a 10x scATAC run are not comparable:

  R1, R3  genomic mates (~50 bp)  -- judged STRICTLY on mean quality and adapter
                                     content; these are what actually align.
  R2      cell barcode (16 bp)    -- judged LENIENTLY: it is short by design, so
                                     length and per-base quality warnings that
                                     matter for a genomic read do not apply. We
                                     only insist it is the expected length.

Deliberately non-fatal, mirroring read_qc.py in the scRNA sibling: chromap's own
adapter handling and the deterministic QC gate downstream are the enforcement
layers. This exists so a human (or the batch reconciler) can see a degraded run
before wondering why its fragment yield was low.
"""
import argparse
import csv
import re
import sys
import zipfile
from pathlib import Path

BARCODE_READ = "R2"          # the 16 bp cell-barcode read
GENOMIC_READS = ("R1", "R3")


def parse_fastqc_zip(path):
    """Pull the few numbers we need out of a FastQC zip's fastqc_data.txt."""
    with zipfile.ZipFile(path) as z:
        name = next((n for n in z.namelist() if n.endswith("fastqc_data.txt")), None)
        if name is None:
            raise ValueError(f"{path} contains no fastqc_data.txt")
        text = z.read(name).decode("utf-8", "replace")

    out = {"total_reads": None, "read_length": None,
           "mean_quality": None, "adapter_fraction": None}

    for line in text.splitlines():
        if line.startswith("Total Sequences"):
            out["total_reads"] = float(line.split("\t")[1])
        elif line.startswith("Sequence length"):
            # "50" or "35-50" -> take the maximum
            val = line.split("\t")[1]
            out["read_length"] = float(val.split("-")[-1])

    # Mean of the per-base mean-quality column.
    block = re.search(r">>Per base sequence quality.*?\n(.*?)>>END_MODULE",
                      text, re.S)
    if block:
        means = []
        for row in block.group(1).splitlines():
            if row.startswith("#") or not row.strip():
                continue
            parts = row.split("\t")
            if len(parts) > 1:
                try:
                    means.append(float(parts[1]))
                except ValueError:
                    pass
        if means:
            out["mean_quality"] = sum(means) / len(means)

    # Peak cumulative adapter content across all adapter columns.
    block = re.search(r">>Adapter Content.*?\n(.*?)>>END_MODULE", text, re.S)
    if block:
        peak = 0.0
        for row in block.group(1).splitlines():
            if row.startswith("#") or not row.strip():
                continue
            parts = row.split("\t")[1:]
            try:
                peak = max(peak, sum(float(p) for p in parts))
            except ValueError:
                pass
        out["adapter_fraction"] = peak / 100.0   # FastQC reports percent

    return out


def judge(read_role, m, args):
    reasons = []
    if m["total_reads"] is not None and m["total_reads"] < args.min_reads:
        reasons.append(f"total_reads={m['total_reads']:.0f}<{args.min_reads}")

    if read_role == BARCODE_READ:
        # Only the length matters for the barcode read.
        if m["read_length"] is not None and m["read_length"] < args.min_barcode_len:
            reasons.append(
                f"barcode_len={m['read_length']:.0f}<{args.min_barcode_len}")
    else:
        if m["mean_quality"] is not None and m["mean_quality"] < args.min_mean_quality:
            reasons.append(
                f"mean_quality={m['mean_quality']:.1f}<{args.min_mean_quality}")
        if (m["adapter_fraction"] is not None
                and m["adapter_fraction"] > args.max_adapter_fraction):
            reasons.append(
                f"adapter={m['adapter_fraction']:.3f}>{args.max_adapter_fraction}")
    return ("PASS" if not reasons else "FLAG"), ";".join(reasons)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fastqc-dir", required=True)
    ap.add_argument("--samples", required=True, help="samples.tsv (for run list)")
    ap.add_argument("--report", required=True)
    ap.add_argument("--passlist", required=True)
    ap.add_argument("--min-reads", type=float, default=1000)
    ap.add_argument("--min-mean-quality", type=float, default=28)
    ap.add_argument("--max-adapter-fraction", type=float, default=0.10)
    ap.add_argument("--min-barcode-len", type=float, default=16)
    args = ap.parse_args()

    zips = sorted(Path(args.fastqc_dir).glob("*_fastqc.zip"))
    rows, passing = [], []
    for z in zips:
        m_role = re.search(r"_(R[123])_fastqc\.zip$", z.name)
        role = m_role.group(1) if m_role else "?"
        srr = z.name.split("_")[0]
        try:
            metrics = parse_fastqc_zip(z)
        except ValueError as e:
            rows.append({"srr": srr, "read": role, "total_reads": None,
                         "read_length": None, "mean_quality": None,
                         "adapter_fraction": None, "status": "FLAG",
                         "reasons": str(e)})
            continue
        status, reason = judge(role, metrics, args)
        rows.append({"srr": srr, "read": role, **metrics,
                     "status": status, "reasons": reason})
        if status == "PASS":
            passing.append(f"{srr}_{role}")

    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        # No FastQC output at all: still write the files so the DAG can proceed.
        Path(args.report).write_text("srr\tread\tstatus\treasons\n")
        Path(args.passlist).write_text("")
        print("read_qc: no FastQC zips found", file=sys.stderr)
        return

    with open(args.report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    Path(args.passlist).write_text("\n".join(passing) + ("\n" if passing else ""))

    n_flag = sum(1 for r in rows if r["status"] == "FLAG")
    print(f"read QC: {len(rows) - n_flag}/{len(rows)} reads clean. "
          f"Report: {args.report}")
    for r in rows:
        if r["status"] == "FLAG":
            print(f"  FLAG {r['srr']} {r['read']}: {r['reasons']}", file=sys.stderr)
    # Non-fatal by design.


if __name__ == "__main__":
    main()

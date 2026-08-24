#!/usr/bin/env python3
"""Completeness reconciler for the scATAC-seq chromap pipeline.

WHY THIS EXISTS
Snakemake treats a sample as "done" when its per-sample .h5ad merely EXISTS.
That is not enough. The object can exist while the sample failed the QC gate,
or while the fragment file it points at has been deleted, or while the peak x
cell matrix is entirely empty. This script enforces the fuller contract and
produces the to-do list of what still needs (re)running.

A sample is DONE iff ALL of:
  1. <workdir>/h5ad/<sample>.h5ad exists and is readable, AND
  2. <sample> is listed in <workdir>/qc/qc_pass.txt (passed the QC gate), AND
  3. uns['fragments'] names a file that still exists (the durable artifact), AND
  4. X has at least one non-zero entry.

Otherwise it is categorized, and each category carries an `action` telling a
driver whether to RERUN it (transient -- hand back to Snakemake, worth retrying)
or FLAG it (genuine -- a rerun of the same inputs cannot fix it, a human decides):

  missing        - no h5ad (not run yet, or an upstream step crashed).    -> RERUN
  corrupt        - h5ad exists but is unreadable (truncated/killed write). -> RERUN
                   (the batch loop quarantines it first so it is rebuilt.)
  no_fragments   - h5ad fine but uns['fragments'] is unset or dangling.
                   Usually keep_fastq=false plus a deleted workdir; the
                   fragments are the durable artifact, so losing them means
                   re-aligning.                                            -> RERUN
  empty_matrix   - h5ad + QC pass but X has zero non-zero entries. The
                   signature of a chrom-naming mismatch between fragments
                   and peaks (chr1 vs 1) -- a wiring fault, not a data one.  -> FLAG
  read_qc_fail   - the raw reads were flagged by read-QC (bad library or
                   wrong chemistry). The same reads re-aligned give the same
                   bad result.                                             -> FLAG
  qc_fail        - h5ad but not in qc_pass.txt (failed the gate); reason
                   carried through from qc_gate.tsv.                       -> FLAG

DIFFERENCES FROM THE scRNA SIBLING
Its `no_layers` category (absent Velocyto layers) has no ATAC analogue -- there
are no spliced/unspliced layers here. In its place: `no_fragments` and
`empty_matrix`, which are the two ways an ATAC object can look complete and be
useless.

EXIT CODES (single-workdir mode)
  0 = every sample DONE
  1 = work remains
  2 = error (bad arguments, unreadable inputs)
Read-only: it launches nothing.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

# category -> what a driver should do about it
ACTIONS = {
    "missing": "RERUN",
    "corrupt": "RERUN",
    "no_fragments": "RERUN",
    "empty_matrix": "FLAG",
    "read_qc_fail": "FLAG",
    "qc_fail": "FLAG",
}


def read_passlist(path):
    p = Path(path)
    return set(p.read_text().split()) if p.exists() else set()


def read_qc_reasons(path):
    """sample -> reason string, from qc_gate.tsv."""
    p = Path(path)
    if not p.exists():
        return {}
    out = {}
    with open(p) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("status") == "FAIL":
                out[row["sample"]] = row.get("reasons", "")
    return out


def flagged_runs(path):
    """runs flagged by read_qc.tsv (non-fatal report layer)."""
    p = Path(path)
    if not p.exists():
        return {}
    out = {}
    with open(p) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("status") == "FLAG":
                out.setdefault(row["srr"], []).append(
                    f"{row.get('read','?')}:{row.get('reasons','')}")
    return out


def inspect(h5ad_path):
    """(ok, problem) for one object. Imports anndata lazily so the reconciler
    still runs (reporting `missing`) in an environment without it."""
    try:
        import anndata as ad
    except ImportError:
        return None, "anndata not importable"
    try:
        a = ad.read_h5ad(h5ad_path)
    except Exception as e:                                   # noqa: BLE001
        return False, f"unreadable: {type(e).__name__}"

    frag = a.uns.get("fragments")
    if not frag or not Path(str(frag)).exists():
        return False, "no_fragments"
    nnz = a.X.nnz if hasattr(a.X, "nnz") else int((a.X != 0).sum())
    if nnz == 0:
        return False, "empty_matrix"
    return True, ""


def reconcile_workdir(workdir, samples_tsv, accession=""):
    workdir = Path(workdir)
    qc = workdir / "qc"
    passing = read_passlist(qc / "qc_pass.txt")
    qc_reasons = read_qc_reasons(qc / "qc_gate.tsv")
    read_flags = flagged_runs(qc / "read_qc.tsv")

    rows = []
    with open(samples_tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    samples = {}
    for r in rows:
        samples.setdefault(r["sample"], []).append(r.get("srr", ""))

    results = []
    for sample, runs in sorted(samples.items()):
        h5 = workdir / "h5ad" / f"{sample}.h5ad"
        rec = {"accession": accession, "sample": sample,
               "h5ad": str(h5), "category": "", "detail": "", "action": ""}

        if not h5.exists():
            rec["category"] = "missing"
            rec["detail"] = "no .h5ad"
        else:
            ok, problem = inspect(h5)
            if ok is None:
                rec["category"] = "missing"
                rec["detail"] = problem
            elif not ok and problem.startswith("unreadable"):
                rec["category"] = "corrupt"
                rec["detail"] = problem
            elif not ok:
                rec["category"] = problem            # no_fragments | empty_matrix
                rec["detail"] = problem
            elif sample not in passing:
                rec["category"] = "qc_fail"
                rec["detail"] = qc_reasons.get(sample, "not in qc_pass.txt")
            else:
                flags = [f"{s}:{','.join(v)}" for s, v in read_flags.items()
                         if s in runs]
                if flags:
                    # DONE, but worth surfacing: the object passed while its raw
                    # reads were flagged. Not a failure -- a caveat.
                    rec["category"] = ""
                    rec["detail"] = "read_qc flags: " + "; ".join(flags)

        rec["action"] = ACTIONS.get(rec["category"], "")
        results.append(rec)

    done = [r for r in results if not r["category"]]
    return {"accession": accession, "workdir": str(workdir),
            "n_samples": len(results), "n_done": len(done),
            "complete": len(done) == len(results) and bool(results),
            "samples": results}


def print_report(res, quiet=False):
    print(f"{res['accession'] or res['workdir']}: "
          f"{res['n_done']}/{res['n_samples']} done")
    if quiet:
        return
    for r in res["samples"]:
        if r["category"]:
            print(f"  {r['action']:5} {r['category']:13} {r['sample']}  {r['detail']}")
        elif r["detail"]:
            print(f"  {'ok':5} {'(caveat)':13} {r['sample']}  {r['detail']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", help="reconcile every study in this manifest TSV")
    ap.add_argument("--samples", help="samples.tsv (single-workdir mode)")
    ap.add_argument("--workdir", help="pipeline workdir (single-workdir mode)")
    ap.add_argument("--accession", default="", help="label for the report")
    ap.add_argument("--report", help="write a per-sample TSV report here")
    ap.add_argument("--json", dest="json_out", help="write machine-readable JSON here")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()

    all_res = []
    if args.manifest:
        with open(args.manifest) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                # An uncommented header row parses as a study named "accession"
                # and produces a phantom entry. Skip it explicitly rather than
                # relying on every manifest being written correctly.
                if f[0].strip().lower() == "accession":
                    continue
                acc, workdir, samples_tsv = f[0], f[1], f[2]
                if not Path(samples_tsv).exists():
                    all_res.append({"accession": acc, "workdir": workdir,
                                    "n_samples": 0, "n_done": 0,
                                    "complete": False, "samples": [
                                        {"accession": acc, "sample": "-",
                                         "category": "missing",
                                         "detail": f"no samples.tsv at {samples_tsv}",
                                         "action": "RERUN", "h5ad": ""}]})
                    continue
                all_res.append(reconcile_workdir(workdir, samples_tsv, acc))
    elif args.samples and args.workdir:
        all_res.append(reconcile_workdir(args.workdir, args.samples, args.accession))
    else:
        ap.error("give --manifest, or both --samples and --workdir")
        sys.exit(2)

    for res in all_res:
        print_report(res, args.quiet)

    if args.report:
        flat = [r for res in all_res for r in res["samples"]]
        Path(args.report).parent.mkdir(parents=True, exist_ok=True)
        with open(args.report, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(flat[0]), delimiter="\t")
            w.writeheader(); w.writerows(flat)
    if args.json_out:
        Path(args.json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json_out).write_text(json.dumps(all_res, indent=2))

    complete = all(r["complete"] for r in all_res) and bool(all_res)
    n_rerun = sum(1 for res in all_res for r in res["samples"]
                  if r["action"] == "RERUN")
    n_flag = sum(1 for res in all_res for r in res["samples"]
                 if r["action"] == "FLAG")
    if not complete:
        print(f"\n{n_rerun} sample(s) to RERUN, {n_flag} needing a human (FLAG)")
    sys.exit(0 if complete else 1)


if __name__ == "__main__":
    main()

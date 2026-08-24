#!/usr/bin/env python3
"""Preflight for the scATAC pipeline: environment, tools, and free disk.

Three checks, all of which must pass before the DAG does any real work:

  1. Python pins  -- the packages the pipeline's Python steps import, compared
     against requirements-pipeline.txt.
  2. Binary tools -- chromap, macs3, fastqc, multiqc on PATH. (prefetch is
     checked but non-fatal: it is only needed for source=sra runs.)
  3. Free disk    -- ATAC FASTQs are large (~30-80 GB per 10x sample) and this
     pipeline is typically run on a shared mount. Aborting here is far cheaper
     than filling a filesystem other people depend on halfway through a study.

Exits non-zero on any hard failure so the workflow stops before downloading.
"""
import argparse
import importlib.metadata as md
import re
import shutil
import sys
from pathlib import Path

# Distribution names the pipeline actually depends on (import name may differ,
# e.g. pyyaml -> yaml). Compared by distribution name against the pins file.
REQUIRED = ["anndata", "numpy", "pandas", "scipy", "h5py",
            "pysam", "snapatac2", "pyranges", "pyyaml"]

# Binaries that must be on PATH. macs3 (not macs2): MACS2 is unmaintained and
# does not build on Python 3.13.
REQUIRED_BINS = ["chromap", "macs3", "fastqc", "multiqc"]
OPTIONAL_BINS = ["prefetch", "fasterq-dump"]   # source=sra runs only

PIN_RE = re.compile(r"^([A-Za-z0-9_.\-]+)==([0-9][^\s\;]*)")


def parse_pins(requirements_path):
    pins = {}
    with open(requirements_path) as fh:
        for line in fh:
            m = PIN_RE.match(line.strip())
            if m:
                pins[m.group(1).lower()] = m.group(2)
    return pins


def installed_version(dist):
    for name in (dist, dist.replace("pyyaml", "PyYAML")):
        try:
            return md.version(name)
        except md.PackageNotFoundError:
            continue
    return None


def check_python(requirements_path):
    pins = parse_pins(requirements_path)
    problems = []
    print(f"{'package':12} {'required':12} {'installed':12} status")
    for pkg in REQUIRED:
        want = pins.get(pkg)
        have = installed_version(pkg)
        if have is None:
            status, _ = "MISSING", problems.append(pkg)
        elif want is None:
            status = "ok (unpinned)"      # present but not pinned; acceptable
        elif have != want:
            status, _ = "MISMATCH", problems.append(pkg)
        else:
            status = "ok"
        print(f"{pkg:12} {str(want):12} {str(have):12} {status}")
    return problems


def check_binaries():
    problems = []
    print(f"\n{'tool':14} {'path':50} status")
    for tool in REQUIRED_BINS:
        path = shutil.which(tool)
        if path is None:
            problems.append(tool)
        print(f"{tool:14} {str(path):50} {'ok' if path else 'MISSING'}")
    for tool in OPTIONAL_BINS:
        path = shutil.which(tool)
        print(f"{tool:14} {str(path):50} "
              f"{'ok' if path else 'absent (only needed for source=sra)'}")
    return problems


def check_disk(check_dir, min_free_gb):
    """Abort if the workdir's filesystem has less than min_free_gb available."""
    if not min_free_gb:
        return []
    target = Path(check_dir)
    # The workdir may not exist yet; walk up to the nearest existing ancestor.
    while not target.exists() and target != target.parent:
        target = target.parent
    usage = shutil.disk_usage(target)
    free_gb = usage.free / 1024 ** 3
    total_gb = usage.total / 1024 ** 3
    ok = free_gb >= min_free_gb
    print(f"\ndisk: {free_gb:,.0f} GB free of {total_gb:,.0f} GB on {target} "
          f"(need {min_free_gb:,} GB) -- {'ok' if ok else 'INSUFFICIENT'}")
    if ok:
        return []
    return [f"only {free_gb:,.0f} GB free on {target}, need {min_free_gb:,} GB"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("requirements", help="path to requirements-pipeline.txt")
    ap.add_argument("--check-dir", default="",
                    help="directory whose filesystem free space is checked")
    ap.add_argument("--min-free-gb", type=float, default=0,
                    help="abort if the check-dir filesystem has less free than this")
    ap.add_argument("--warn-only", action="store_true",
                    help="report problems but exit 0")
    args = ap.parse_args()

    problems = check_python(args.requirements)
    problems += check_binaries()
    if args.check_dir:
        problems += check_disk(args.check_dir, args.min_free_gb)

    if problems and not args.warn_only:
        print(f"\nPreflight FAILED: {problems}", file=sys.stderr)
        print("Fix the above, then rerun. See README 'Install'.", file=sys.stderr)
        sys.exit(1)
    print("\nPreflight OK" if not problems else "\n(warn-only) problems ignored")


if __name__ == "__main__":
    main()

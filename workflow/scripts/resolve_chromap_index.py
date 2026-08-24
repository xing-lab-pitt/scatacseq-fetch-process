#!/usr/bin/env python3
"""Produce the chromap index the aligner will use -- rebuilding on mismatch.

ONE PRINCIPLE: the index and the aligner must come from the same environment.
The chromap rule can only read the pipeline-owned path this script writes, so a
version mismatch between "the index someone built last year" and "the chromap on
PATH today" is impossible by construction.

Unlike STAR, chromap does not stamp a readable format version into its index, so
we keep a sidecar `<index>.build_info` recording the chromap version that built
it. Decision table:

  supplied path exists + build_info matches running chromap  -> symlink it (instant)
  supplied path exists + build_info missing or mismatched    -> rebuild into it
  supplied path empty                                        -> build ephemerally
                                                                in the workdir

This mirrors resolve_star_index.py in the scRNA sibling, which probes STAR's
genomeVersion for the same purpose.
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path


def chromap_version():
    """chromap prints its version to stderr on --version; tolerate either stream."""
    try:
        p = subprocess.run(["chromap", "--version"],
                           capture_output=True, text=True, check=False)
    except FileNotFoundError:
        print("chromap not found on PATH. See README 'Install' -- chromap is not "
              "a Python package and must be built from source.", file=sys.stderr)
        sys.exit(1)
    return (p.stdout.strip() or p.stderr.strip() or "unknown").splitlines()[0]


def build_info_path(index_path):
    return Path(str(index_path) + ".build_info")


def build_index(fasta, dest, version):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Building chromap index: {fasta} -> {dest}  (chromap {version})")
    subprocess.run(["chromap", "-i", "-r", str(fasta), "-o", str(dest)], check=True)
    build_info_path(dest).write_text(f"{version}\n{fasta}\n")
    print(f"Built {dest} ({dest.stat().st_size:,} bytes)")


def link(src, out):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.is_symlink() or out.exists():
        out.unlink()
    out.symlink_to(Path(src).resolve())
    print(f"Symlinked {out} -> {Path(src).resolve()}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True,
                    help="pipeline-owned index path the aligner will read")
    ap.add_argument("--supplied", default="",
                    help="prebuilt index hint, or a writable path to build+keep at")
    ap.add_argument("--fasta", required=True, help="reference FASTA")
    args = ap.parse_args()

    version = chromap_version()
    fasta = Path(args.fasta)
    if not fasta.exists():
        print(f"Reference FASTA not found: {fasta}", file=sys.stderr)
        sys.exit(1)

    supplied = Path(args.supplied) if args.supplied else None

    if supplied is None:
        # No location given: build ephemerally at the pipeline-owned path.
        build_index(fasta, args.out, version)
        return

    if supplied.exists() and supplied.stat().st_size > 0:
        info = build_info_path(supplied)
        recorded = info.read_text().splitlines()[0].strip() if info.exists() else None
        if recorded == version:
            link(supplied, args.out)
            return
        print(f"Supplied index {supplied} was built by "
              f"{recorded or 'an unrecorded chromap version'}, but the chromap on "
              f"PATH is {version}. Rebuilding to keep index and aligner in step.")

    build_index(fasta, supplied, version)
    link(supplied, args.out)


if __name__ == "__main__":
    main()

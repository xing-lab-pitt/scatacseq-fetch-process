"""The orientation detector is the pipeline's cheapest safety net: it turns a
6-CPU-hour silent failure (an almost-empty fragment file) into a 6-second one.
These tests pin the three outcomes that matter."""
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import revcomp, write_fastq

SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/detect_barcode_orientation.py"


def run(fastq, whitelist, out_wl, report, extra=()):
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--fastq", str(fastq),
         "--whitelist", str(whitelist), "--out-whitelist", str(out_wl),
         "--report", str(report), *extra],
        capture_output=True, text=True)


def parse_report(path):
    return dict(line.split("\t", 1)
                for line in Path(path).read_text().splitlines() if "\t" in line)


def test_forward_reads_choose_forward(tmp_path, whitelist, rng):
    wl_path, codes = whitelist
    fq = tmp_path / "bc.fastq.gz"
    write_fastq(fq, [rng.choice(codes) for _ in range(2000)])

    r = run(fq, wl_path, tmp_path / "out.txt", tmp_path / "rep.txt")
    assert r.returncode == 0, r.stderr

    rep = parse_report(tmp_path / "rep.txt")
    assert rep["orientation_chosen"] == "forward"
    assert float(rep["forward_match_rate"]) > 0.95
    # The emitted whitelist must be byte-identical in content to the input.
    assert set((tmp_path / "out.txt").read_text().split()) == set(codes)


def test_revcomp_reads_choose_revcomp(tmp_path, whitelist, rng):
    wl_path, codes = whitelist
    fq = tmp_path / "bc.fastq.gz"
    write_fastq(fq, [revcomp(rng.choice(codes)) for _ in range(2000)])

    r = run(fq, wl_path, tmp_path / "out.txt", tmp_path / "rep.txt")
    assert r.returncode == 0, r.stderr

    rep = parse_report(tmp_path / "rep.txt")
    assert rep["orientation_chosen"] == "revcomp"
    assert float(rep["revcomp_match_rate"]) > 0.95
    # Every emitted barcode must map back into the original whitelist.
    emitted = (tmp_path / "out.txt").read_text().split()
    assert all(revcomp(b) in set(codes) for b in emitted)


def test_unmatched_barcodes_abort(tmp_path, whitelist, rng):
    """Random barcodes match neither orientation -> hard stop, not a quiet pass.

    This is the case that matters most: without it chromap runs to completion and
    emits a near-empty fragment file that looks like a successful run."""
    wl_path, _ = whitelist
    fq = tmp_path / "bc.fastq.gz"
    write_fastq(fq, ["".join(rng.choice("ACGT") for _ in range(16))
                     for _ in range(2000)])

    r = run(fq, wl_path, tmp_path / "out.txt", tmp_path / "rep.txt")
    assert r.returncode == 1
    assert "Neither orientation matches" in r.stderr
    # The report is still written, so a human can see both rates.
    assert (tmp_path / "rep.txt").exists()


def test_forced_orientation_overrides_detection(tmp_path, whitelist, rng):
    """--orientation revcomp must win even when forward matches perfectly."""
    wl_path, codes = whitelist
    fq = tmp_path / "bc.fastq.gz"
    write_fastq(fq, [rng.choice(codes) for _ in range(2000)])

    r = run(fq, wl_path, tmp_path / "out.txt", tmp_path / "rep.txt",
            extra=("--orientation", "revcomp"))
    assert r.returncode == 0, r.stderr
    rep = parse_report(tmp_path / "rep.txt")
    assert rep["orientation_chosen"] == "revcomp"
    assert rep["orientation_mode"] == "revcomp"

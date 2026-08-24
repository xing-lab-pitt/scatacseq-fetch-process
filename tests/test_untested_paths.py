"""Paths that existed but had never been exercised.

Each test here corresponds to a line in VERIFICATION_BACKLOG.md. They are grouped
in one file because what they share is provenance, not subject: every one covers
code that was written, reviewed, believed correct, and never actually run.
"""
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow/scripts"))

from fragments_to_h5ad import nucleosome_signal, MONO_NUCLEOSOME  # noqa: E402
from prepare_runs import md5s_for  # noqa: E402

CHECK_VERSIONS = ROOT / "workflow/scripts/check_versions.py"


# --------------------------- nucleosome signal -------------------------------
# Computed by hand rather than by snapATAC2 (whose frag_size_distr is aggregate
# only), so it had no coverage at all -- it appeared in check_h5ad's required
# column list and nowhere else.

def write_fragments(path, rows):
    """rows: (chrom, start, end, barcode)"""
    with gzip.open(path, "wt") as fh:
        for c, s, e, bc in rows:
            fh.write(f"{c}\t{s}\t{e}\t{bc}\t1\n")


def test_nucleosome_signal_is_long_over_short(tmp_path):
    """Ratio of mono-nucleosomal (>=147 bp) to sub-nucleosomal fragments."""
    f = tmp_path / "frag.tsv.gz"
    rows = []
    rows += [("chr1", i * 1000, i * 1000 + 80, "BC1") for i in range(10)]    # short
    rows += [("chr1", i * 1000, i * 1000 + 200, "BC1") for i in range(30)]   # long
    write_fragments(f, rows)
    got = nucleosome_signal(f, ["BC1"])
    assert got[0] == pytest.approx(30 / 10)


def test_nucleosome_signal_boundary_is_inclusive(tmp_path):
    """A fragment of exactly MONO_NUCLEOSOME counts as long, not short."""
    f = tmp_path / "frag.tsv.gz"
    write_fragments(f, [("chr1", 0, MONO_NUCLEOSOME, "BC1"),
                        ("chr1", 100, 100 + MONO_NUCLEOSOME - 1, "BC1")])
    got = nucleosome_signal(f, ["BC1"])
    assert got[0] == pytest.approx(1.0)   # one long, one short


def test_nucleosome_signal_is_per_barcode(tmp_path):
    """Barcodes must not contaminate each other."""
    f = tmp_path / "frag.tsv.gz"
    rows = [("chr1", i * 500, i * 500 + 200, "LONG") for i in range(20)]
    rows += [("chr1", i * 500, i * 500 + 80, "SHORT") for i in range(20)]
    write_fragments(f, rows)
    got = nucleosome_signal(f, ["LONG", "SHORT"])
    assert got[0] > 10          # all long / floor of 1
    assert got[1] == pytest.approx(0.0)


def test_nucleosome_signal_ignores_unrequested_barcodes(tmp_path):
    """Memory is bounded by the number of CELLS, not the whitelist, so barcodes
    outside the requested set must be skipped entirely."""
    f = tmp_path / "frag.tsv.gz"
    write_fragments(f, [("chr1", 0, 200, "WANTED"), ("chr1", 0, 200, "OTHER")])
    got = nucleosome_signal(f, ["WANTED"])
    assert len(got) == 1


def test_nucleosome_signal_no_short_fragments_is_finite(tmp_path):
    """A barcode with zero sub-nucleosomal fragments must not divide by zero."""
    f = tmp_path / "frag.tsv.gz"
    write_fragments(f, [("chr1", i * 500, i * 500 + 300, "BC1") for i in range(5)])
    got = nucleosome_signal(f, ["BC1"])
    assert got[0] == pytest.approx(5.0)     # 5 long / floor of 1
    assert got[0] == got[0]                 # not NaN


# ------------------------------ disk guard -----------------------------------
# min_free_gb exists to stop a multi-sample study filling a shared mount at 99%
# capacity. It had never actually aborted anything.

def test_disk_guard_aborts_when_space_is_short(tmp_path):
    req = ROOT / "requirements-pipeline.txt"
    r = subprocess.run(
        [sys.executable, str(CHECK_VERSIONS), str(req),
         "--check-dir", str(tmp_path), "--min-free-gb", "99999999"],
        capture_output=True, text=True)
    assert r.returncode == 1
    assert "INSUFFICIENT" in r.stdout
    assert "free on" in r.stderr or "Preflight FAILED" in r.stderr


def test_disk_guard_passes_when_space_is_ample(tmp_path):
    req = ROOT / "requirements-pipeline.txt"
    r = subprocess.run(
        [sys.executable, str(CHECK_VERSIONS), str(req),
         "--check-dir", str(tmp_path), "--min-free-gb", "0.001"],
        capture_output=True, text=True)
    assert "-- ok" in r.stdout


def test_disk_guard_walks_up_to_an_existing_parent(tmp_path):
    """The workdir usually does not exist yet at preflight time, so the check
    must resolve to the nearest existing ancestor rather than erroring."""
    missing = tmp_path / "not" / "created" / "yet"
    r = subprocess.run(
        [sys.executable, str(CHECK_VERSIONS), str(ROOT / "requirements-pipeline.txt"),
         "--check-dir", str(missing), "--min-free-gb", "0.001"],
        capture_output=True, text=True)
    assert "disk:" in r.stdout


# ------------------------------ md5 selection --------------------------------
# The md5 columns were written empty on every row, so download_fastq's checksum
# verification could never fire.

def test_md5s_align_with_chosen_urls():
    urls = ["a_1.fq.gz", "a_2.fq.gz", "a_3.fq.gz"]
    md5_map = {"SRR1": ["md5_one", "md5_two", "md5_three"]}
    # chosen order is (genomic1, barcode, genomic2) = (_1, _2, _3) here
    assert md5s_for(md5_map, "SRR1", urls, (urls[0], urls[1], urls[2])) == \
        ("md5_one", "md5_two", "md5_three")


def test_md5s_follow_reordered_urls():
    """Roles are assigned by READ LENGTH, so the barcode may not be the middle
    file. The checksums must follow the files, not the slot order."""
    urls = ["a_1.fq.gz", "a_2.fq.gz", "a_3.fq.gz"]
    md5_map = {"SRR1": ["m1", "m2", "m3"]}
    # barcode turned out to be _1; genomic reads are _2 and _3
    assert md5s_for(md5_map, "SRR1", urls, (urls[1], urls[0], urls[2])) == \
        ("m2", "m1", "m3")


def test_md5s_absent_yields_empty_not_wrong():
    """ENA does not always publish checksums. Better an empty column (check
    skipped) than a mismatched one (every download fails)."""
    assert md5s_for({}, "SRR1", ["a", "b", "c"], ("a", "b", "c")) == ("", "", "")


def test_md5s_count_mismatch_yields_empty():
    """If ENA lists a different number of checksums than files, the alignment is
    unknowable -- refuse to guess."""
    md5_map = {"SRR1": ["m1", "m2"]}          # 2 checksums, 3 files
    assert md5s_for(md5_map, "SRR1", ["a", "b", "c"], ("a", "b", "c")) == ("", "", "")

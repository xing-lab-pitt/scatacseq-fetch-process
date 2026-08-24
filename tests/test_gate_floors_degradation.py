"""Gate-floor validation: prove the QC floors CATCH bad data, not just pass good.

WHY THIS FILE EXISTS
Every sample run through this pipeline so far has been healthy, which only
demonstrates the floors do not produce false FAILs. It says nothing about whether
they would catch a broken sample -- a floor of FRiP >= 0.15 is useless if broken
data still scores 0.4. A threshold that has only ever been shown data that passes
is "known permissive", not "known discriminating".

So: take the metric profile of a REAL passing sample (atac_pbmc_500_nextgem,
measured), degrade one axis at a time to values characteristic of a known failure
mode, and assert the gate fails for the RIGHT reason. Each case names the real
laboratory or bioinformatic failure it stands for.

These are metric-level tests. The upstream input-level guards (barcode
orientation, read geometry, contig naming) are covered in their own test files --
together they are defence in depth, and the orientation check in particular fires
long before the gate would ever see the sample.
"""
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/qc_gate.py"

# Measured on the real reference run: atac_pbmc_500_nextgem, 10% downsample,
# GRCh38, blacklist-filtered peaks. This is a genuinely good sample.
REAL_GOOD = dict(
    frac_reads_in_peaks=0.574,
    median_tss_enrichment=21.4,
    n_cells=451,
    valid_barcode_frac=0.939,
    duplicate_rate=0.170,
)

# The production floors, from config.example.yaml.
FLOORS = dict(
    min_frac_reads_in_peaks=0.15,
    min_tss_enrichment=4.0,
    min_estimated_cells=100,
    min_valid_barcode_frac=0.70,
    max_duplicate_rate=0.80,
)


def make_h5ad(path, **overrides):
    metrics = {**REAL_GOOD, **overrides}
    a = ad.AnnData(
        X=sparse.csr_matrix(np.ones((4, 3))),
        obs=pd.DataFrame(index=[f"bc{i}" for i in range(4)]),
        var=pd.DataFrame(index=[f"chr1:{i}00-{i}50" for i in range(3)]),
    )
    a.uns.update(metrics)
    a.write_h5ad(path)
    return path


def run_gate(tmp_path, h5ad):
    args = [sys.executable, str(SCRIPT), f"S={h5ad}",
            "--report", str(tmp_path / "r.tsv"),
            "--passlist", str(tmp_path / "p.txt")]
    for k, v in FLOORS.items():
        args += [f"--{k.replace('_', '-')}", str(v)]
    r = subprocess.run(args, capture_output=True, text=True)
    return r, (tmp_path / "p.txt").read_text().split()


def test_the_real_good_sample_passes(tmp_path):
    """Baseline. If this fails the floors are too strict and everything below
    is meaningless."""
    _, passed = run_gate(tmp_path, make_h5ad(tmp_path / "s.h5ad"))
    assert passed == ["S"]


# --- one axis at a time, each standing for a real failure mode ---------------

@pytest.mark.parametrize("override,expect_metric,failure_mode", [
    # Tn5 never inserted efficiently, or peaks called on too few fragments:
    # signal is spread across the genome rather than concentrated in peaks.
    (dict(frac_reads_in_peaks=0.05), "frac_reads_in_peaks",
     "failed transposition / dead chromatin"),
    # Degraded or over-fixed nuclei: no enrichment at promoters, which is the
    # single most diagnostic ATAC quality signal.
    (dict(median_tss_enrichment=1.2), "median_tss_enrichment",
     "degraded nuclei, no promoter enrichment"),
    # Loading or barcode failure: almost nothing passes cell calling.
    (dict(n_cells=12), "n_cells",
     "loading failure / wrong whitelist"),
    # Wrong whitelist or wrong orientation: most barcodes are not real.
    (dict(valid_barcode_frac=0.31), "valid_barcode_frac",
     "wrong whitelist or wrong barcode read"),
    # Over-amplified low-complexity library: nearly all reads are duplicates.
    (dict(duplicate_rate=0.94), "duplicate_rate",
     "over-amplified low-complexity library"),
])
def test_single_axis_degradation_is_caught(tmp_path, override, expect_metric,
                                           failure_mode):
    r, passed = run_gate(tmp_path, make_h5ad(tmp_path / "s.h5ad", **override))
    assert passed == [], f"gate PASSED a sample representing: {failure_mode}"
    assert expect_metric in r.stderr, (
        f"gate failed but blamed the wrong metric for: {failure_mode}\n{r.stderr}")


def test_catastrophic_sample_fails_on_several_axes(tmp_path):
    """A truly broken run degrades on many axes at once; the report should name
    them all, so a human sees the pattern rather than one symptom."""
    r, passed = run_gate(tmp_path, make_h5ad(
        tmp_path / "s.h5ad",
        frac_reads_in_peaks=0.02, median_tss_enrichment=0.9,
        n_cells=3, valid_barcode_frac=0.05, duplicate_rate=0.99))
    assert passed == []
    for metric in ("frac_reads_in_peaks", "median_tss_enrichment", "n_cells",
                   "valid_barcode_frac", "duplicate_rate"):
        assert metric in r.stderr


# --- where the floors actually sit, relative to real data --------------------

@pytest.mark.parametrize("metric,floor_key,good,direction", [
    ("frac_reads_in_peaks", "min_frac_reads_in_peaks", 0.574, "min"),
    ("median_tss_enrichment", "min_tss_enrichment", 21.4, "min"),
    ("n_cells", "min_estimated_cells", 451, "min"),
])
def test_floor_is_below_real_data_but_not_absurdly_so(metric, floor_key, good,
                                                      direction):
    """Documents the margin between the floor and real good data.

    This is not a pass/fail quality bar -- it is a regression guard on a
    JUDGEMENT. The floors sit 3-5x below observed good values, which is
    deliberate for a floor, but if someone later raises one to within 20% of
    real data it will start rejecting healthy samples, and this test says so.
    """
    floor = FLOORS[floor_key]
    ratio = good / floor
    assert ratio > 1.5, (
        f"{metric}: floor {floor} is within 1.5x of real good data ({good}) -- "
        "too close, healthy samples will start failing")
    assert ratio < 20, (
        f"{metric}: floor {floor} is {ratio:.0f}x below real good data ({good}) "
        "-- so permissive it may not catch anything; revisit")


def test_just_below_floor_fails_and_just_above_passes(tmp_path):
    """The floor discriminates at its own boundary, in both directions."""
    _, passed_below = run_gate(
        tmp_path, make_h5ad(tmp_path / "lo.h5ad", frac_reads_in_peaks=0.149))
    assert passed_below == []
    _, passed_above = run_gate(
        tmp_path, make_h5ad(tmp_path / "hi.h5ad", frac_reads_in_peaks=0.151))
    assert passed_above == ["S"]

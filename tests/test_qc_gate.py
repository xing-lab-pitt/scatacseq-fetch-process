"""The QC gate must be deterministic and its thresholds must come from config,
never from the data. These tests pin the boundary behaviour, which is where a
gate silently inverts."""
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/qc_gate.py"


def make_h5ad(path, *, frip=0.5, tsse=10.0, n_cells=500,
              valid_bc=0.9, dup=0.3):
    n_obs, n_var = 10, 5
    a = ad.AnnData(
        X=sparse.csr_matrix(np.ones((n_obs, n_var))),
        obs=pd.DataFrame(index=[f"bc{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=[f"chr1:{i}00-{i}50" for i in range(n_var)]),
    )
    a.uns.update({
        "frac_reads_in_peaks": frip,
        "median_tss_enrichment": tsse,
        "n_cells": n_cells,
        "valid_barcode_frac": valid_bc,
        "duplicate_rate": dup,
    })
    a.write_h5ad(path)
    return path


def run_gate(tmp_path, entries, **thresholds):
    args = [sys.executable, str(SCRIPT), *entries,
            "--report", str(tmp_path / "r.tsv"),
            "--passlist", str(tmp_path / "p.txt")]
    for k, v in thresholds.items():
        args += [f"--{k.replace('_', '-')}", str(v)]
    r = subprocess.run(args, capture_output=True, text=True)
    passed = (tmp_path / "p.txt").read_text().split()
    return r, passed


DEFAULTS = dict(min_frac_reads_in_peaks=0.15, min_tss_enrichment=4.0,
                min_estimated_cells=100, min_valid_barcode_frac=0.70,
                max_duplicate_rate=0.80)


def test_good_sample_passes(tmp_path):
    make_h5ad(tmp_path / "S.h5ad")
    r, passed = run_gate(tmp_path, [f"S={tmp_path/'S.h5ad'}"], **DEFAULTS)
    assert passed == ["S"]


def test_low_frip_fails_and_names_the_metric(tmp_path):
    make_h5ad(tmp_path / "S.h5ad", frip=0.05)
    r, passed = run_gate(tmp_path, [f"S={tmp_path/'S.h5ad'}"], **DEFAULTS)
    assert passed == []
    assert "frac_reads_in_peaks" in r.stderr


def test_high_duplicate_rate_fails(tmp_path):
    """duplicate_rate is a MAX threshold -- the one direction that inverts easily."""
    make_h5ad(tmp_path / "S.h5ad", dup=0.95)
    r, passed = run_gate(tmp_path, [f"S={tmp_path/'S.h5ad'}"], **DEFAULTS)
    assert passed == []
    assert "duplicate_rate" in r.stderr


def test_value_exactly_at_threshold_passes(tmp_path):
    """>= for minimums, <= for maximums. Pinned so a later refactor can't flip it."""
    make_h5ad(tmp_path / "S.h5ad", frip=0.15, dup=0.80)
    r, passed = run_gate(tmp_path, [f"S={tmp_path/'S.h5ad'}"], **DEFAULTS)
    assert passed == ["S"]


def test_missing_metric_fails_rather_than_passes(tmp_path):
    """A metric that was never computed must NOT be treated as acceptable."""
    p = tmp_path / "S.h5ad"
    make_h5ad(p)
    a = ad.read_h5ad(p)
    del a.uns["median_tss_enrichment"]
    a.write_h5ad(p)
    r, passed = run_gate(tmp_path, [f"S={p}"], **DEFAULTS)
    assert passed == []
    assert "median_tss_enrichment:missing" in r.stderr


def test_gate_always_exits_zero(tmp_path):
    """Reporting is separate from enforcement: the report must survive failure."""
    make_h5ad(tmp_path / "S.h5ad", frip=0.0)
    r, _ = run_gate(tmp_path, [f"S={tmp_path/'S.h5ad'}"], **DEFAULTS)
    assert r.returncode == 0
    assert (tmp_path / "r.tsv").exists()


def test_mixed_batch_partitions_correctly(tmp_path):
    make_h5ad(tmp_path / "good.h5ad")
    make_h5ad(tmp_path / "bad.h5ad", tsse=1.0)
    r, passed = run_gate(
        tmp_path,
        [f"good={tmp_path/'good.h5ad'}", f"bad={tmp_path/'bad.h5ad'}"],
        **DEFAULTS)
    assert passed == ["good"]

"""check_h5ad answers 'is this object well-formed?' before the QC gate asks
'is this good data?'. A malformed object with plausible numbers is the failure
mode that costs the most downstream time."""
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/check_h5ad.py"

OBS_COLS = ["n_fragments", "frip", "tss_enrichment",
            "nucleosome_signal", "frac_mito", "is_cell"]


def make(path, fragments_file, *, n_obs=5, n_var=4, empty_x=False,
         drop_obs=None, nan_col=None, bad_frag=False):
    X = sparse.csr_matrix((n_obs, n_var)) if empty_x else \
        sparse.csr_matrix(np.ones((n_obs, n_var)))
    obs = pd.DataFrame({c: np.ones(n_obs) for c in OBS_COLS},
                       index=[f"bc{i}" for i in range(n_obs)])
    if drop_obs:
        obs = obs.drop(columns=[drop_obs])
    if nan_col:
        obs[nan_col] = np.nan
    var = pd.DataFrame(
        {"chrom": ["chr1"] * n_var,
         "start": list(range(n_var)),
         "end": [i + 50 for i in range(n_var)]},
        index=[f"chr1:{i}00-{i}50" for i in range(n_var)])
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.uns["fragments"] = "/nonexistent/frags.tsv.gz" if bad_frag \
        else str(fragments_file)
    a.write_h5ad(path)
    return path


def run(path, tmp_path):
    return subprocess.run(
        [sys.executable, str(SCRIPT), str(path), "--out", str(tmp_path / "ok")],
        capture_output=True, text=True)


def test_well_formed_object_passes(tmp_path):
    frag = tmp_path / "f.tsv.gz"; frag.write_bytes(b"")
    r = run(make(tmp_path / "a.h5ad", frag), tmp_path)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "ok").exists()


def test_empty_matrix_is_caught(tmp_path):
    """All-zero X is the classic chrom-name-mismatch symptom."""
    frag = tmp_path / "f.tsv.gz"; frag.write_bytes(b"")
    r = run(make(tmp_path / "a.h5ad", frag, empty_x=True), tmp_path)
    assert r.returncode == 1
    assert "zero non-zero entries" in r.stderr


def test_missing_obs_column_is_caught(tmp_path):
    frag = tmp_path / "f.tsv.gz"; frag.write_bytes(b"")
    r = run(make(tmp_path / "a.h5ad", frag, drop_obs="frip"), tmp_path)
    assert r.returncode == 1
    assert "frip" in r.stderr


def test_all_nan_metric_is_caught(tmp_path):
    """A column that was allocated but never computed must not pass."""
    frag = tmp_path / "f.tsv.gz"; frag.write_bytes(b"")
    r = run(make(tmp_path / "a.h5ad", frag, nan_col="tss_enrichment"), tmp_path)
    assert r.returncode == 1
    assert "entirely NaN" in r.stderr


def test_dangling_fragments_path_is_caught(tmp_path):
    frag = tmp_path / "f.tsv.gz"; frag.write_bytes(b"")
    r = run(make(tmp_path / "a.h5ad", frag, bad_frag=True), tmp_path)
    assert r.returncode == 1
    assert "does not exist" in r.stderr


def test_failure_does_not_write_marker(tmp_path):
    frag = tmp_path / "f.tsv.gz"; frag.write_bytes(b"")
    run(make(tmp_path / "a.h5ad", frag, empty_x=True), tmp_path)
    assert not (tmp_path / "ok").exists()

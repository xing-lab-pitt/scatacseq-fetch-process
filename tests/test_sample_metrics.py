"""Sample-level QC metrics must describe the CELLS, not the ambient background.

WHY THIS FILE EXISTS
The gate tests (test_gate_floors_degradation.py) inject pre-made metric values
and check the thresholds behave. They passed while the metrics themselves were
computed wrongly -- like testing that a fire alarm triggers at 100 C without ever
checking the thermometer. This file tests the thermometer.

THE BUG IT PINS
Sample metrics were originally computed over ALL barcodes. In droplet data most
barcodes are empty, so a median over all of them measures the background. Worse,
the error GROWS WITH SEQUENCING DEPTH: deeper runs push more empty droplets above
the import floor, dragging the median further into the noise. Measured on two
real datasets:

                     cells/barcodes   median TSSe (all)   median TSSe (cells)
  atac_pbmc_500         451/594  76%        21.45               22.74
  atac_v1_pbmc_5k     5026/63614   8%         1.36               23.37

Both are good samples. On the all-barcode median the deeper one scores 17x worse
and would FAIL a floor of 4.0. The metric was tracking depth, not quality.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from fragments_to_h5ad import sample_metrics  # noqa: E402

CELL_TSSE, AMBIENT_TSSE = 23.0, 1.2
CELL_FRIP, AMBIENT_FRIP = 0.75, 0.06


def make_obs(n_cells, n_ambient):
    """A realistic droplet mixture: a few real cells among many empty droplets."""
    rows = []
    for _ in range(n_cells):
        rows.append(dict(n_fragments=8000.0, frip=CELL_FRIP,
                         tss_enrichment=CELL_TSSE, frac_dup=0.30, is_cell=True))
    for _ in range(n_ambient):
        rows.append(dict(n_fragments=150.0, frip=AMBIENT_FRIP,
                         tss_enrichment=AMBIENT_TSSE, frac_dup=0.05, is_cell=False))
    return pd.DataFrame(rows)


def test_metrics_describe_cells_not_ambient():
    m = sample_metrics(make_obs(500, 5000))
    assert m["median_tss_enrichment"] == pytest.approx(CELL_TSSE)
    assert m["frac_reads_in_peaks"] == pytest.approx(CELL_FRIP)
    assert m["n_cells"] == 500


def test_metric_is_invariant_to_ambient_count():
    """THE regression guard.

    Adding empty droplets must not move the score. This is exactly what broke:
    the same library sequenced deeper gained ambient barcodes and its TSSe
    'dropped' from 21 to 1.4 without the cells changing at all.
    """
    shallow = sample_metrics(make_obs(500, 200))     # 71% cells, like PBMC500
    deep = sample_metrics(make_obs(500, 50_000))     # 1% cells, deeper run

    assert shallow["median_tss_enrichment"] == pytest.approx(
        deep["median_tss_enrichment"]), \
        "TSS enrichment changed when only the ambient droplet count changed"
    assert shallow["frac_reads_in_peaks"] == pytest.approx(
        deep["frac_reads_in_peaks"])
    assert shallow["n_cells"] == deep["n_cells"] == 500


def test_the_old_all_barcode_median_would_have_failed():
    """Demonstrates the bug is real, so the guard above cannot be dismissed."""
    obs = make_obs(500, 50_000)
    all_barcode_median = float(np.nanmedian(obs["tss_enrichment"]))
    cells_only = sample_metrics(obs)["median_tss_enrichment"]

    assert all_barcode_median == pytest.approx(AMBIENT_TSSE)   # ~1.2
    assert cells_only == pytest.approx(CELL_TSSE)              # ~23
    # A floor of 4.0 rejects the good sample under the old computation.
    assert all_barcode_median < 4.0 < cells_only


def test_frip_is_depth_weighted_within_cells():
    """A cell with 50x the fragments should count 50x, so shallow cells cannot
    drag a good sample down."""
    obs = pd.DataFrame([
        dict(n_fragments=50_000.0, frip=0.80, tss_enrichment=20.0,
             frac_dup=0.3, is_cell=True),
        dict(n_fragments=1_000.0, frip=0.20, tss_enrichment=20.0,
             frac_dup=0.3, is_cell=True),
    ])
    m = sample_metrics(obs)
    unweighted = (0.80 + 0.20) / 2
    weighted = (0.80 * 50_000 + 0.20 * 1_000) / 51_000
    assert m["frac_reads_in_peaks"] == pytest.approx(weighted)
    assert m["frac_reads_in_peaks"] != pytest.approx(unweighted)


def test_ambient_figures_are_retained_for_diagnosis():
    """The background numbers are useful -- for spotting a FAILED cell call --
    they just must not be what the gate reads."""
    m = sample_metrics(make_obs(100, 1000))
    assert m["ambient_median_tss_enrichment"] == pytest.approx(AMBIENT_TSSE)
    assert m["median_tss_enrichment"] == pytest.approx(CELL_TSSE)
    assert m["cell_fraction_of_barcodes"] == pytest.approx(100 / 1100)


def test_no_cells_called_yields_no_metrics_rather_than_zeros():
    """If cell calling produced nothing, the metrics must be absent (which the
    gate treats as a FAIL) rather than silently ambient values that might pass."""
    m = sample_metrics(make_obs(0, 5000))
    assert m["n_cells"] == 0
    assert m["median_tss_enrichment"] is None
    assert m["duplicate_rate"] is None
    assert m["frac_reads_in_peaks"] == 0.0


def test_missing_tsse_is_none_not_nan():
    """have_tsse=False must yield None, so the gate reports 'missing' and fails,
    rather than a NaN that could slip through a comparison."""
    m = sample_metrics(make_obs(50, 50), have_tsse=False)
    assert m["median_tss_enrichment"] is None

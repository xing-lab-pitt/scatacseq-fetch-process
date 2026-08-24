"""Multiome (ARC) machinery: translation table + MuData pairing.

These run entirely on synthetic whitelists and synthetic objects, so all of the
multiome logic is testable BEFORE any real multiome FASTQ exists. That matters
because the real failure mode is silent: a wrong translation table still yields
a complete, well-formed, QC-passing object keyed to the wrong cells. The tests
below pin the guards that turn that silence into an error.
"""
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
TRANSLATE = ROOT / "workflow/scripts/resolve_arc_translation.py"
MERGE = ROOT / "workflow/scripts/merge_h5ad.py"


def write_barcodes(path, codes):
    path.write_text("\n".join(codes) + "\n")
    return path


def make_pair(tmp_path, rng, n=500):
    """A synthetic ARC whitelist pair: equal length, disjoint, both 16 bp."""
    def draw(prefix):
        seen, out = set(), []
        while len(out) < n:
            bc = prefix + "".join(rng.choice("ACGT") for _ in range(14))
            if bc not in seen:
                seen.add(bc); out.append(bc)
        return out
    atac = write_barcodes(tmp_path / "atac.txt", draw("AA"))
    gex = write_barcodes(tmp_path / "gex.txt", draw("GG"))
    return atac, gex


def run_translate(atac, gex, out):
    return subprocess.run(
        [sys.executable, str(TRANSLATE), "--atac-whitelist", str(atac),
         "--gex-whitelist", str(gex), "-o", str(out)],
        capture_output=True, text=True)


# --------------------------- translation table -------------------------------

def test_translation_table_pairs_by_line(tmp_path, rng):
    atac, gex = make_pair(tmp_path, rng)
    out = tmp_path / "t.tsv"
    r = run_translate(atac, gex, out)
    assert r.returncode == 0, r.stderr

    a_lines = atac.read_text().split()
    g_lines = gex.read_text().split()
    table = [l.split("\t") for l in out.read_text().splitlines()]
    assert len(table) == len(a_lines)
    # Row N pairs line N of each list -- but written TARGET FIRST, i.e.
    # GEX<TAB>ATAC. chromap keys on column 2 and emits column 1
    # (src/barcode_translator.h). This assertion previously encoded the
    # opposite, wrong order; the integration test caught it.
    for i, (to_bc, from_bc) in enumerate(table):
        assert from_bc == a_lines[i], "column 2 must be the ATAC (source) barcode"
        assert to_bc == g_lines[i], "column 1 must be the GEX (target) barcode"


def test_unequal_line_counts_refused(tmp_path, rng):
    """Pairing is BY LINE NUMBER, so unequal files mean every mapping is wrong."""
    atac, gex = make_pair(tmp_path, rng)
    short = tmp_path / "short.txt"
    short.write_text("\n".join(gex.read_text().split()[:-5]) + "\n")
    r = run_translate(atac, short, tmp_path / "t.tsv")
    assert r.returncode == 1
    assert "line counts differ" in r.stderr


def test_same_file_twice_refused(tmp_path, rng):
    """The realistic mistake: both ARC whitelists are named 737K-arc-v1.txt."""
    atac, _ = make_pair(tmp_path, rng)
    r = run_translate(atac, atac, tmp_path / "t.tsv")
    assert r.returncode == 1
    assert "IDENTICAL" in r.stderr


def test_overlapping_universes_refused(tmp_path, rng):
    """Real ARC lists are disjoint; heavy overlap means a file is not what it
    claims to be."""
    atac, gex = make_pair(tmp_path, rng)
    codes = atac.read_text().split()
    mixed = codes[:400] + gex.read_text().split()[400:]
    write_barcodes(tmp_path / "mixed.txt", mixed)
    r = run_translate(atac, tmp_path / "mixed.txt", tmp_path / "t.tsv")
    assert r.returncode == 1
    assert "BOTH lists" in r.stderr


def test_wrong_barcode_length_refused(tmp_path, rng):
    atac, gex = make_pair(tmp_path, rng)
    write_barcodes(tmp_path / "short_bc.txt",
                   [c[:12] for c in gex.read_text().split()])
    r = run_translate(atac, tmp_path / "short_bc.txt", tmp_path / "t.tsv")
    assert r.returncode == 1
    assert "16 bp" in r.stderr


# ----------------------------- MuData pairing --------------------------------

def make_obj(path, barcodes, n_var=5, var_prefix="chr1:"):
    a = ad.AnnData(
        X=sparse.csr_matrix(np.ones((len(barcodes), n_var))),
        obs=pd.DataFrame(index=list(barcodes)),
        var=pd.DataFrame(index=[f"{var_prefix}{i}00-{i}50" for i in range(n_var)]),
    )
    a.uns["fragments"] = str(path) + ".frag"
    a.write_h5ad(path)
    return path


def run_merge(atac_h5ad, out, gex_h5ad=None):
    args = [sys.executable, str(MERGE), "--h5ad", str(atac_h5ad), "-o", str(out)]
    if gex_h5ad:
        args += ["--gex-h5ad", str(gex_h5ad)]
    return subprocess.run(args, capture_output=True, text=True)


def test_translated_barcodes_pair_into_mudata(tmp_path, rng):
    """After translation the ATAC object is keyed by GEX barcodes, so the two
    objects share obs_names and pair cleanly."""
    pytest.importorskip("mudata")
    shared = [f"GG{''.join(rng.choice('ACGT') for _ in range(14))}" for _ in range(200)]
    atac = make_obj(tmp_path / "atac.h5ad", shared)
    gex = make_obj(tmp_path / "gex.h5ad", shared, var_prefix="GENE")
    r = run_merge(atac, tmp_path / "out.h5mu", gex)
    assert r.returncode == 0, r.stderr
    assert "100.0%" in r.stdout
    assert (tmp_path / "out.h5mu").exists()


def test_untranslated_barcodes_refuse_to_pair(tmp_path, rng):
    """THE test for this feature. Without --barcode-translate the two objects are
    keyed to different barcode universes. Both are individually well-formed and
    would each pass their own QC gate -- only the intersection reveals it."""
    pytest.importorskip("mudata")
    atac_bcs = [f"AA{''.join(rng.choice('ACGT') for _ in range(14))}" for _ in range(200)]
    gex_bcs = [f"GG{''.join(rng.choice('ACGT') for _ in range(14))}" for _ in range(200)]
    atac = make_obj(tmp_path / "atac.h5ad", atac_bcs)
    gex = make_obj(tmp_path / "gex.h5ad", gex_bcs, var_prefix="GENE")

    r = run_merge(atac, tmp_path / "out.h5mu", gex)
    assert r.returncode == 1
    assert "FAILED pairing" in r.stderr
    assert "barcode-translate" in r.stderr
    # And it must NOT leave a plausible-looking artifact behind.
    assert not (tmp_path / "out.h5mu").exists()


def test_partial_overlap_below_threshold_refused(tmp_path, rng):
    pytest.importorskip("mudata")
    shared = [f"GG{''.join(rng.choice('ACGT') for _ in range(14))}" for _ in range(40)]
    only_atac = [f"AA{''.join(rng.choice('ACGT') for _ in range(14))}" for _ in range(160)]
    only_gex = [f"CC{''.join(rng.choice('ACGT') for _ in range(14))}" for _ in range(160)]
    atac = make_obj(tmp_path / "atac.h5ad", shared + only_atac)
    gex = make_obj(tmp_path / "gex.h5ad", shared + only_gex, var_prefix="GENE")
    r = run_merge(atac, tmp_path / "out.h5mu", gex)
    assert r.returncode == 1          # 40/200 = 20%, below the 30% floor


def test_atac_only_merge_still_writes_h5ad(tmp_path, rng):
    """The plain 10x path must be untouched by any of the multiome machinery."""
    bcs = [f"AA{''.join(rng.choice('ACGT') for _ in range(14))}" for _ in range(50)]
    atac = make_obj(tmp_path / "atac.h5ad", bcs)
    out = tmp_path / "combined.h5ad"
    r = run_merge(atac, out)
    assert r.returncode == 0, r.stderr
    assert out.exists()
    assert ad.read_h5ad(out).n_obs == 50

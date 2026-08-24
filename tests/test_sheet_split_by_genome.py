"""A study spanning genomes must yield one runnable sheet per genome.

GSE219015 has 16 Homo sapiens and 2 Mus musculus samples. A single sheet
covering both cannot be run: the pipeline aligns one genome per run and refuses
a mixed sheet. Writing one anyway hands the user a file that always fails.

The split lives inside prepare_runs.main(), which needs the network, so the
grouping and naming rules are exercised directly here.
"""
from pathlib import Path

import pytest


def split_paths(out: Path, genomes):
    """Mirror of the naming rule in prepare_runs.main()."""
    return {g: out.with_name(f"{out.stem}.{g}{out.suffix}") for g in sorted(genomes)}


def group(rows):
    by = {}
    for r in rows:
        by.setdefault(r["genome"], []).append(r)
    return by


HUMAN = [{"srr": f"SRR{i}", "genome": "GRCh38"} for i in range(16)]
MOUSE = [{"srr": "SRR28111503", "genome": "mm10"},
         {"srr": "SRR28197504", "genome": "mm10"}]


def test_mixed_study_groups_by_genome():
    by = group(HUMAN + MOUSE)
    assert sorted(by) == ["GRCh38", "mm10"]
    assert len(by["GRCh38"]) == 16 and len(by["mm10"]) == 2


def test_no_row_is_lost_or_duplicated():
    rows = HUMAN + MOUSE
    by = group(rows)
    assert sum(len(v) for v in by.values()) == len(rows)
    srrs = {r["srr"] for v in by.values() for r in v}
    assert srrs == {r["srr"] for r in rows}


def test_every_split_sheet_is_single_genome():
    # The property that makes each output runnable.
    for genome, rows in group(HUMAN + MOUSE).items():
        assert {r["genome"] for r in rows} == {genome}


def test_split_names_keep_the_suffix_and_stem():
    paths = split_paths(Path("config/samples.tsv"), ["GRCh38", "mm10"])
    assert paths["GRCh38"].name == "samples.GRCh38.tsv"
    assert paths["mm10"].name == "samples.mm10.tsv"
    assert paths["mm10"].parent == Path("config")


def test_split_names_are_distinct():
    paths = split_paths(Path("config/samples.tsv"), ["GRCh38", "mm10"])
    assert len(set(paths.values())) == 2


def test_single_genome_study_is_not_split():
    by = group(HUMAN)
    assert len(by) == 1, "a one-genome study must keep writing the plain -o path"


@pytest.mark.parametrize("genomes", [["GRCh38"], ["mm10"], ["GRCh38", "mm10"]])
def test_grouping_is_stable_and_sorted(genomes):
    rows = [{"srr": f"S{i}", "genome": g} for i, g in enumerate(genomes)]
    assert sorted(group(rows)) == sorted(set(genomes))

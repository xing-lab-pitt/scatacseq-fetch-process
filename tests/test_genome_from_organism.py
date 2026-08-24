"""Genome must follow the organism SRA reports, not the --genome default.

GSE219015 is Mus musculus. The sheet recorded GRCh38 because prepare_runs.py
never read the species. 97% of reads failed to map, and the run still produced a
complete object with FRiP 0.70 and TSSe 21 that passed the QC gate, because
every downstream metric is computed from the reads that did map.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
from prepare_runs import ORGANISM_GENOME, genome_for  # noqa: E402


def run(organism, srr="SRR28197504"):
    return {"srr": srr, "organism": organism}


def test_mouse_gives_mm10():
    assert genome_for(run("Mus musculus"), "GRCh38", explicit=False) == "mm10"


def test_human_gives_grch38():
    assert genome_for(run("Homo sapiens"), "GRCh38", explicit=False) == "GRCh38"


def test_organism_beats_the_default():
    # The real failure: default GRCh38, mouse data, no --genome passed.
    assert genome_for(run("Mus musculus"), "GRCh38", explicit=False) != "GRCh38"


def test_case_and_whitespace_tolerated():
    assert genome_for(run("  mus musculus "), "GRCh38", explicit=False) == "mm10"


def test_explicit_genome_contradicting_organism_is_refused():
    with pytest.raises(SystemExit) as e:
        genome_for(run("Mus musculus"), "GRCh38", explicit=True)
    msg = str(e.value)
    assert "Mus musculus" in msg and "mm10" in msg


def test_explicit_genome_agreeing_with_organism_is_fine():
    assert genome_for(run("Mus musculus"), "mm10", explicit=True) == "mm10"


def test_unknown_organism_falls_back_to_default():
    assert genome_for(run("Danio rerio"), "GRCh38", explicit=False) == "GRCh38"


def test_missing_organism_falls_back_to_default():
    assert genome_for({"srr": "SRR1"}, "GRCh38", explicit=False) == "GRCh38"


def test_every_mapped_genome_is_one_the_pipeline_can_resolve():
    # A genome string that no config references block knows would fail at parse
    # time in the Snakefile, which is a worse place to find out.
    assert set(ORGANISM_GENOME.values()) <= {"GRCh38", "mm10"}

"""prepare_runs must say "unknown" rather than guess "10x".

The modality decides which whitelist and barcode offset are used. It used to
default to "10x" whenever the metadata was silent, which reads downstream
exactly like a confirmed value. GSE219015's mouse sample is multiome and no
metadata field says so, so it was labelled 10x with no sign anything was
uncertain.

"unknown" is honest: detect_barcode_orientation then measures it from the reads.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
from prepare_runs import (  # noqa: E402
    MULTIOME_BARCODE_READ_LEN,
    barcode_len_in_prose,
    detect_modality,
)

GEO_MOUSE = ("atac library (cell barcodes are matchable with mrna, mtdna, "
             "targetsite, and cellhash lib in sample 15): read1 file contains "
             "one end of accessible region; i5 file (1-24) contains cell "
             "barcodes, read2 file contains the other end")
GEO_HUMAN = ("libraries were prepared following 10x genomics chromium next gem "
             "single cell multiome atac + gene expression user manual")
PLAIN_10X = ("libraries were prepared with the chromium single cell atac kit; "
             "i5 file (1-16) contains cell barcodes")


def run(**kw):
    base = {"srr": "SRR1", "gsm": "GSM1", "library_name": "GSM1",
            "sample_name": "GSM1", "library_strategy": "ATAC-seq"}
    base.update(kw)
    return base


# --- barcode length parsed out of prose -------------------------------------

def test_reads_24bp_barcode_from_geo_wording():
    assert barcode_len_in_prose(GEO_MOUSE) == MULTIOME_BARCODE_READ_LEN


def test_reads_16bp_barcode_from_geo_wording():
    assert barcode_len_in_prose(PLAIN_10X) == 16


def test_no_barcode_length_stated():
    assert barcode_len_in_prose("no geometry mentioned anywhere here") is None


def test_unrelated_parenthetical_is_not_a_barcode_length():
    assert barcode_len_in_prose("samples (1-24) were collected in batches") is None


# --- modality from prose -----------------------------------------------------

def test_named_kit_gives_multiome(monkeypatch):
    monkeypatch.setattr("prepare_runs.library_prose", lambda r: GEO_HUMAN)
    assert detect_modality(run()) == "multiome"


def test_24bp_geometry_gives_multiome_without_the_word(monkeypatch):
    # The mouse sample: never says "multiome", but states a 24 bp barcode.
    monkeypatch.setattr("prepare_runs.library_prose", lambda r: GEO_MOUSE)
    assert detect_modality(run()) == "multiome"


def test_16bp_geometry_gives_10x(monkeypatch):
    monkeypatch.setattr("prepare_runs.library_prose", lambda r: PLAIN_10X)
    assert detect_modality(run()) == "10x"


def test_silent_metadata_gives_unknown_not_10x(monkeypatch):
    monkeypatch.setattr("prepare_runs.library_prose",
                        lambda r: "bone marrow was collected from donors")
    assert detect_modality(run()) == "unknown"


def test_empty_prose_gives_unknown(monkeypatch):
    monkeypatch.setattr("prepare_runs.library_prose", lambda r: "")
    assert detect_modality(run()) == "unknown"


def test_offline_gives_unknown_not_10x():
    assert detect_modality(run(), use_network=False) == "unknown"


def test_name_field_still_wins_without_network():
    assert detect_modality(run(library_name="sample_ARC-v1"),
                           use_network=False) == "multiome"


@pytest.mark.parametrize("phrase", [
    "cell barcodes are matchable with mrna lib",
    "joint single-cell rna and atac data",
])
def test_rna_partner_implies_multiome(monkeypatch, phrase):
    monkeypatch.setattr("prepare_runs.library_prose", lambda r: phrase)
    assert detect_modality(run()) == "multiome"


def test_10x_is_never_returned_on_a_guess(monkeypatch):
    # The whole point: "10x" must be evidence-backed, never a fallback.
    monkeypatch.setattr("prepare_runs.library_prose", lambda r: "nothing useful")
    assert detect_modality(run()) != "10x"

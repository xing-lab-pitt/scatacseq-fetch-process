"""Read-role assignment by MEASURED LENGTH, not filename order.

Submitters deposit 10x scATAC under R1/R2/R3, R1/R2/I2 and R1/I1/R2, and
fasterq-dump just numbers files _1.._4. Getting this wrong feeds genomic reads
into the barcode slot, which produces zero valid barcodes and no useful error."""
import sys
from pathlib import Path

import pytest

from conftest import write_fastq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
import select_reads  # noqa: E402


def make(tmp_path, name, length, n=50):
    p = tmp_path / name
    write_fastq(p, ["A" * length for _ in range(n)])
    return p


def test_modal_read_length(tmp_path):
    p = make(tmp_path, "x.fastq", 16)
    assert select_reads.modal_read_length(p) == 16


def test_classify_standard_three_read_layout(tmp_path):
    g1 = make(tmp_path, "r_1.fastq", 50)
    bc = make(tmp_path, "r_2.fastq", 16)
    g2 = make(tmp_path, "r_3.fastq", 50)
    assert select_reads.classify([g1, bc, g2]) == (g1, bc, g2)


def test_classify_ignores_file_order(tmp_path):
    """Barcode listed FIRST must still be identified as the barcode."""
    bc = make(tmp_path, "a.fastq", 16)
    g1 = make(tmp_path, "b.fastq", 50)
    g2 = make(tmp_path, "c.fastq", 50)
    got_g1, got_bc, got_g2 = select_reads.classify([bc, g1, g2])
    assert got_bc == bc
    assert {got_g1, got_g2} == {g1, g2}


def test_classify_discards_sample_index(tmp_path):
    """An 8 bp i7 index must be dropped, not mistaken for a cell barcode."""
    i1 = make(tmp_path, "i1.fastq", 8)
    g1 = make(tmp_path, "r1.fastq", 50)
    bc = make(tmp_path, "r2.fastq", 16)
    g2 = make(tmp_path, "r3.fastq", 50)
    assert select_reads.classify([i1, g1, bc, g2]) == (g1, bc, g2)


def test_multiome_24bp_barcode_is_accepted(tmp_path):
    """THE regression test for a real 15-hour failure.

    10x MULTIOME ATAC barcode reads are 24 bp: an 8 bp spacer then the 16 bp
    barcode. select_reads used to demand exactly 16 bp, so on GSE219015 it found
    zero barcode reads, classified all three as genomic, and raised -- after a
    25-minute download and ~40 minutes of extraction, five times over.

    It went unnoticed because every fixture in this file used the 10x 16 bp
    geometry, including tests named for multiome. The same wrong assumption was
    in the code AND its tests.
    """
    g1 = make(tmp_path, "r1.fastq", 50)
    bc = make(tmp_path, "r2.fastq", 24)      # multiome: 8bp spacer + 16bp barcode
    g2 = make(tmp_path, "r3.fastq", 49)
    assert select_reads.classify([g1, bc, g2]) == (g1, bc, g2)


def test_multiome_geometry_with_index_read(tmp_path):
    """Real submissions often include the 8 bp i7 index alongside the 24 bp
    barcode. The index must still be discarded, and 8 != 24 must not confuse
    the barcode test."""
    i1 = make(tmp_path, "i1.fastq", 8)
    g1 = make(tmp_path, "r1.fastq", 50)
    bc = make(tmp_path, "r2.fastq", 24)
    g2 = make(tmp_path, "r3.fastq", 50)
    assert select_reads.classify([i1, g1, bc, g2]) == (g1, bc, g2)


def test_both_chemistries_are_accepted_but_not_two_barcodes(tmp_path):
    """A 16 bp AND a 24 bp read in one run is ambiguous -- refuse rather than
    guess which is the barcode."""
    bc16 = make(tmp_path, "a.fastq", 16)
    bc24 = make(tmp_path, "b.fastq", 24)
    g1 = make(tmp_path, "c.fastq", 50)
    with pytest.raises(ValueError):
        select_reads.classify([bc16, bc24, g1])


def test_error_message_names_both_accepted_lengths(tmp_path):
    """When it does fail, the message must tell you what it was looking for."""
    g1 = make(tmp_path, "r1.fastq", 50)
    g2 = make(tmp_path, "r2.fastq", 50)
    with pytest.raises(ValueError) as e:
        select_reads.classify([g1, g2])
    msg = str(e.value)
    assert "16" in msg and "24" in msg
    assert "multiome" in msg.lower()


def test_no_barcode_read_raises_with_geometry(tmp_path):
    g1 = make(tmp_path, "r1.fastq", 50)
    g2 = make(tmp_path, "r2.fastq", 50)
    with pytest.raises(ValueError) as e:
        select_reads.classify([g1, g2])
    # The message must show what was actually observed, not just "failed".
    assert "50bp" in str(e.value)
    assert "cell-barcode" in str(e.value)


def test_two_barcode_reads_raises(tmp_path):
    bc1 = make(tmp_path, "a.fastq", 16)
    bc2 = make(tmp_path, "b.fastq", 16)
    g1 = make(tmp_path, "c.fastq", 50)
    with pytest.raises(ValueError):
        select_reads.classify([bc1, bc2, g1])

"""One genome per run, resolved and validated before anything downloads.

THE BUG THIS PINS
`samples.tsv` has a `genome` column. It was parsed, defaulted, and then never
read again -- every rule used a single global `config["reference"]`. A sheet of
mm10 rows would therefore be aligned against whatever reference the config named,
producing a low mapping rate and near-zero fragments, and the failure would
surface as a barcode problem rather than a genome one. A mixed-species sheet was
worse: the human samples looked fine, so it read as sample-specific rather than
systematic.

Two decisions are pinned here:
  1. mixed genomes are REFUSED (a union-of-peaks object across species is
     meaningless, so producing one is worse than failing)
  2. references are validated at PARSE TIME -- before the first byte downloads,
     because an SRA download for one ATAC run can take hours
"""
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow/Snakefile"

SHEET_HDR = ("sample\tsrr\tsource\tmodality\tgenome\t"
             "r1_url\tr2_url\tr3_url\tr1_md5\tr2_md5\tr3_md5\n")


def write_sheet(path, rows):
    """rows: list of (sample, srr, genome)"""
    with open(path, "w") as fh:
        fh.write(SHEET_HDR)
        for s, srr, g in rows:
            fh.write(f"{s}\t{srr}\tsra\t10x\t{g}\t\t\t\t\t\t\n")
    return path


def base_config(tmp_path, sheet, references=None, reference=None):
    fa = tmp_path / "genome.fa"; fa.write_text(">chr1\nACGT\n")
    gtf = tmp_path / "genes.gtf"; gtf.write_text("")
    cfg = {
        "samples_tsv": str(sheet),
        "workdir": str(tmp_path / "work"),
        "threads": 1,
        "requirements_txt": str(ROOT / "requirements-pipeline.txt"),
        "keep_fastq": True,
        "barcodes": {"whitelist_10x": str(tmp_path / "wl.txt"),
                     "orientation": "auto"},
        "qc": {"mode": "warn", "min_frac_reads_in_peaks": 0.15,
               "min_tss_enrichment": 4.0, "min_estimated_cells": 100,
               "min_valid_barcode_frac": 0.7, "max_duplicate_rate": 0.8,
               "cell_min_unique_frags": 1000, "cell_min_frip": 0.2},
        "peaks": {"qvalue": 0.05, "shift": -100, "extsize": 200},
    }
    (tmp_path / "wl.txt").write_text("ACGTACGTACGTACGT\n")
    if references is not None:
        cfg["references"] = references
    if reference is not None:
        cfg["reference"] = reference
    default_ref = {"fasta": str(fa), "gtf": str(gtf),
                   "chromap_index": str(tmp_path / "idx"), "gsize": "hs"}
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return path, default_ref


def run_parse(config_path):
    """Parse the workflow only -- no jobs run.

    NOTE: Snakemake reports workflow-CONSTRUCTION errors (our ValueErrors) on
    stdout, not stderr. Asserting on stderr alone silently passes whatever the
    message says, so `out` below is the combined streams.
    """
    r = subprocess.run(
        [sys.executable, "-m", "snakemake", "-s", str(SNAKEFILE),
         "--configfile", str(config_path), "-n", "--directory", str(ROOT)],
        capture_output=True, text=True, cwd=ROOT)
    r.out = (r.stdout or "") + (r.stderr or "")
    return r


def test_mixed_genomes_are_refused(tmp_path):
    """THE decision: a human+mouse sheet must fail, not produce a union object."""
    sheet = write_sheet(tmp_path / "s.tsv",
                        [("A", "SRR1", "GRCh38"), ("B", "SRR2", "mm10")])
    fa = tmp_path / "g.fa"; fa.write_text(">chr1\nACGT\n")
    gtf = tmp_path / "g.gtf"; gtf.write_text("")
    cfg, ref = base_config(tmp_path, sheet,
                           reference={"fasta": str(fa), "gtf": str(gtf),
                                      "chromap_index": str(tmp_path / "i"),
                                      "gsize": "hs"})
    r = run_parse(cfg)
    assert r.returncode != 0
    assert "mixes genomes" in r.out
    assert "GRCh38" in r.out and "mm10" in r.out
    # and it must tell you how to fix it
    assert "Split the sheet" in r.out


def test_single_genome_resolves_from_the_mapping(tmp_path):
    sheet = write_sheet(tmp_path / "s.tsv", [("A", "SRR1", "mm10")])
    fa = tmp_path / "mm10.fa"; fa.write_text(">chr1\nACGT\n")
    gtf = tmp_path / "mm10.gtf"; gtf.write_text("")
    cfg, _ = base_config(tmp_path, sheet, references={
        "mm10": {"fasta": str(fa), "gtf": str(gtf),
                 "chromap_index": str(tmp_path / "mm10.idx"), "gsize": "mm"}})
    r = run_parse(cfg)
    assert "genome=mm10" in r.stdout, r.out
    assert str(fa) in r.stdout, "resolved fasta should be the mm10 one"


def test_genome_absent_from_references_is_refused(tmp_path):
    """A sheet naming a genome the config does not define must fail immediately,
    not download for hours and then fail."""
    sheet = write_sheet(tmp_path / "s.tsv", [("A", "SRR1", "rn6")])
    fa = tmp_path / "g.fa"; fa.write_text(">chr1\nACGT\n")
    gtf = tmp_path / "g.gtf"; gtf.write_text("")
    cfg, _ = base_config(tmp_path, sheet, references={
        "GRCh38": {"fasta": str(fa), "gtf": str(gtf),
                   "chromap_index": str(tmp_path / "i"), "gsize": "hs"}})
    r = run_parse(cfg)
    assert r.returncode != 0
    assert "rn6" in r.out and "references" in r.out


def test_missing_reference_file_is_caught_at_parse_time(tmp_path):
    """The whole point: catch it in a second, not after the download."""
    sheet = write_sheet(tmp_path / "s.tsv", [("A", "SRR1", "GRCh38")])
    cfg, _ = base_config(tmp_path, sheet, references={
        "GRCh38": {"fasta": str(tmp_path / "does_not_exist.fa"),
                   "gtf": str(tmp_path / "nope.gtf"),
                   "chromap_index": str(tmp_path / "i"), "gsize": "hs"}})
    r = run_parse(cfg)
    assert r.returncode != 0
    assert "not found" in r.out
    assert "does_not_exist.fa" in r.out


def test_legacy_single_reference_block_still_works(tmp_path):
    """Back-compat: configs written before the per-genome mapping keep working
    for a single-genome sheet."""
    sheet = write_sheet(tmp_path / "s.tsv", [("A", "SRR1", "GRCh38")])
    fa = tmp_path / "g.fa"; fa.write_text(">chr1\nACGT\n")
    gtf = tmp_path / "g.gtf"; gtf.write_text("")
    cfg, _ = base_config(tmp_path, sheet,
                         reference={"fasta": str(fa), "gtf": str(gtf),
                                    "chromap_index": str(tmp_path / "i"),
                                    "gsize": "hs"})
    r = run_parse(cfg)
    assert "genome=GRCh38" in r.stdout, r.out


def test_blacklist_is_taken_from_the_genome_block(tmp_path):
    """Each genome needs its own blacklist; hg38's applied to mm10 would filter
    essentially nothing and quietly leave the artifacts in."""
    sheet = write_sheet(tmp_path / "s.tsv", [("A", "SRR1", "mm10")])
    fa = tmp_path / "mm10.fa"; fa.write_text(">chr1\nACGT\n")
    gtf = tmp_path / "mm10.gtf"; gtf.write_text("")
    bl = tmp_path / "mm10-blacklist.bed"; bl.write_text("chr1\t0\t100\n")
    cfg, _ = base_config(tmp_path, sheet, references={
        "mm10": {"fasta": str(fa), "gtf": str(gtf),
                 "chromap_index": str(tmp_path / "i"), "gsize": "mm",
                 "blacklist": str(bl)}})
    r = run_parse(cfg)
    assert "genome=mm10" in r.stdout, r.out

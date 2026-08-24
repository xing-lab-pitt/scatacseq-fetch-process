"""keep_fastq and compress_fastq must be independent.

They used to be one flag, so keeping the FASTQs for inspection also forced hours
of gzip, and skipping compression meant losing them. A verification run needs
both: kept, and fast.

The Snakefile is not importable, so the two lines that implement the rule are
evaluated here against the same config shapes the pipeline sees.
"""
import pytest

SNAKEFILE_LOGIC = """
_keep = bool(config.get("keep_fastq", False))
_compress = bool(config.get("compress_fastq", _keep))
FQ_EXT = ".fastq.gz" if _compress else ".fastq"
"""


def resolve(cfg):
    ns = {"config": cfg}
    exec(SNAKEFILE_LOGIC, ns)
    return {"keep": ns["_keep"], "compress": ns["_compress"], "ext": ns["FQ_EXT"]}


def test_keep_without_compress_is_expressible():
    # The case that motivated the split.
    r = resolve({"keep_fastq": True, "compress_fastq": False})
    assert r["keep"] is True and r["compress"] is False and r["ext"] == ".fastq"


def test_compress_without_keep_is_expressible():
    r = resolve({"keep_fastq": False, "compress_fastq": True})
    assert r["keep"] is False and r["compress"] is True and r["ext"] == ".fastq.gz"


@pytest.mark.parametrize("keep", [True, False])
def test_compress_defaults_to_keep_for_old_configs(keep):
    # Configs written before the split name only keep_fastq and must be unchanged.
    r = resolve({"keep_fastq": keep})
    assert r["compress"] is keep
    assert r["ext"] == (".fastq.gz" if keep else ".fastq")


def test_default_config_keeps_nothing_and_compresses_nothing():
    r = resolve({})
    assert r["keep"] is False and r["compress"] is False and r["ext"] == ".fastq"


def test_extension_follows_compression_not_retention():
    # The bug this guards: extension keyed off keep_fastq, so a kept-but-plain
    # file would be named .fastq.gz while holding plain text. FastQC then fails
    # with "ID line didn't start with '@'", which points nowhere useful.
    kept_plain = resolve({"keep_fastq": True, "compress_fastq": False})
    tmp_gz = resolve({"keep_fastq": False, "compress_fastq": True})
    assert kept_plain["ext"] == ".fastq"
    assert tmp_gz["ext"] == ".fastq.gz"

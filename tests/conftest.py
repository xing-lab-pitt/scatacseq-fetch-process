"""Shared synthetic fixtures. No network, no cluster, no real data."""
import gzip
import random
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def write_fastq(path, seqs):
    opener = gzip.open if str(path).endswith(".gz") else open
    mode = "wt" if str(path).endswith(".gz") else "w"
    with opener(path, mode) as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@read{i}\n{s}\n+\n{'I' * len(s)}\n")


@pytest.fixture
def rng():
    return random.Random(1234)


@pytest.fixture
def whitelist(tmp_path, rng):
    """A small barcode whitelist file plus the list of barcodes in it."""
    codes = ["".join(rng.choice("ACGT") for _ in range(16)) for _ in range(2000)]
    path = tmp_path / "whitelist.txt"
    path.write_text("\n".join(codes) + "\n")
    return path, codes

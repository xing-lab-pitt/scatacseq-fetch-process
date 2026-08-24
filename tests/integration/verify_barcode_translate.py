#!/usr/bin/env python3
"""INTEGRATION TEST: does chromap's --barcode-translate do what we assume?

THE UNCERTAINTY THIS RESOLVES
chromap documents the flag as only:
    --barcode-translate FILE   Convert barcode to the specified sequences during output
It does not state the file format, the column order, or the direction. The
pipeline originally assumed `source<TAB>target` (ATAC first). THAT GUESS WAS
WRONG, and this test caught it: chromap aborts with "Barcode does not exist in
the translation table." The real format, from src/barcode_translator.h, is
`target<TAB>source` -- GEX first, ATAC second; chromap keys on column 2 and
emits column 1.

Why the guess is dangerous: a wrong format most likely does NOT crash chromap.
It emits untranslated (ATAC) barcodes, so the ATAC and GEX objects share no
barcodes and cannot pair -- and that looks like a failed experiment rather than
a wrong file format.

THE DESIGN
Synthetic, so ground truth is exact:
  - a tiny genome we control
  - reads carrying barcodes taken from the REAL ARC ATAC whitelist, at known
    line numbers
  - the REAL translation table built by resolve_arc_translation.py
Then: run chromap with --barcode-translate and check what comes out.

  PASS  -> output barcodes are the GEX barcodes at the SAME line numbers
  FAIL  -> output barcodes are still the ATAC ones (flag silently ignored, or
           our column order is backwards)

Requires: chromap on PATH, and the two real ARC whitelists.
"""
import gzip
import random
import subprocess
import sys
from pathlib import Path

import os

# Point this at the directory holding the 10x ARC whitelists. They are ~12 MB
# each and licensed by 10x, so they are not in the repo:
#     export SCATAC_WHITELIST_DIR=/path/to/10x_whitelists
WL = Path(os.environ.get("SCATAC_WHITELIST_DIR", "reference/10x_whitelists"))
ATAC_WL = WL / "737K-arc-v1.ATAC.txt"
GEX_WL = WL / "737K-arc-v1.txt"
TABLE = WL / "737K-arc-v1.translation.tsv"

if not ATAC_WL.exists():
    print(f"SKIP: whitelists not found under {WL}. "
          "Set SCATAC_WHITELIST_DIR to the directory holding 737K-arc-v1*.txt.")
    raise SystemExit(0)

N_CELLS = 120
READS_PER_CELL = 300
READ_LEN = 50
GENOME_LEN = 1_200_000


def build_genome(path, rng):
    seq = "".join(rng.choice("ACGT") for _ in range(GENOME_LEN))
    with open(path, "w") as fh:
        fh.write(">chr1\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    return seq


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def main():
    work = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/bctranslate_test")
    work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20260819)

    for f in (ATAC_WL, GEX_WL, TABLE):
        if not f.exists():
            print(f"MISSING required input: {f}", file=sys.stderr)
            sys.exit(2)

    # --- ground truth: the first N_CELLS lines of each whitelist ------------
    atac_bcs = ATAC_WL.read_text().split()[:N_CELLS]
    gex_bcs = GEX_WL.read_text().split()[:N_CELLS]
    expected = dict(zip(atac_bcs, gex_bcs))
    print(f"ground truth: {len(expected)} ATAC->GEX pairs from whitelist lines 1-{N_CELLS}")
    print(f"  e.g. {atac_bcs[0]} -> {gex_bcs[0]}")

    # --- synthetic genome + reads ------------------------------------------
    genome = build_genome(work / "genome.fa", rng)
    r1, r2, r3 = [], [], []
    for bc in atac_bcs:
        for _ in range(READS_PER_CELL):
            start = rng.randint(0, GENOME_LEN - 400)
            flen = rng.choice([120, 200, 260])
            frag = genome[start:start + flen]
            r1.append(frag[:READ_LEN])
            r3.append(revcomp(frag[-READ_LEN:]))
            r2.append(bc)                      # ATAC barcode goes IN

    def write_fq(path, seqs):
        with gzip.open(path, "wt") as fh:
            for i, s in enumerate(seqs):
                fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")

    write_fq(work / "R1.fastq.gz", r1)
    write_fq(work / "R2.fastq.gz", r2)
    write_fq(work / "R3.fastq.gz", r3)
    print(f"synthetic reads: {len(r1):,} triples across {N_CELLS} cells")

    # --- chromap index ------------------------------------------------------
    idx = work / "genome.index"
    subprocess.run(["chromap", "-i", "-r", str(work / "genome.fa"), "-o", str(idx)],
                   check=True, capture_output=True)

    def run_chromap(out_bed, translate):
        cmd = ["chromap", "-t", "4", "--preset", "atac",
               "-x", str(idx), "-r", str(work / "genome.fa"),
               "-1", str(work / "R1.fastq.gz"),
               "-2", str(work / "R3.fastq.gz"),
               "-b", str(work / "R2.fastq.gz"),
               "--barcode-whitelist", str(ATAC_WL),
               "-o", str(out_bed)]
        if translate:
            cmd += ["--barcode-translate", str(TABLE)]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            print(f"chromap failed (translate={bool(translate)}):\n{p.stderr[-1500:]}",
                  file=sys.stderr)
            sys.exit(1)
        return {l.split("\t")[3] for l in
                Path(out_bed).read_text().splitlines() if l.strip()}

    # --- control: no translation -> expect ATAC barcodes out ----------------
    print("\n--- control run: WITHOUT --barcode-translate ---")
    plain = run_chromap(work / "plain.bed", translate=False)
    plain_atac = len(plain & set(atac_bcs))
    plain_gex = len(plain & set(gex_bcs))
    print(f"  distinct output barcodes: {len(plain)}")
    print(f"    matching ATAC whitelist: {plain_atac}")
    print(f"    matching GEX  whitelist: {plain_gex}")

    # --- the actual test: with translation -> expect GEX barcodes out -------
    print("\n--- test run: WITH --barcode-translate ---")
    translated = run_chromap(work / "translated.bed", translate=True)
    tr_atac = len(translated & set(atac_bcs))
    tr_gex = len(translated & set(gex_bcs))
    print(f"  distinct output barcodes: {len(translated)}")
    print(f"    matching ATAC whitelist: {tr_atac}")
    print(f"    matching GEX  whitelist: {tr_gex}")

    # --- verdict ------------------------------------------------------------
    print("\n" + "=" * 68)
    if not translated:
        print("VERDICT: FAIL - translated run produced no fragments at all.")
        sys.exit(1)

    if tr_gex > tr_atac and tr_gex >= 0.8 * len(expected):
        # And confirm the mapping is the RIGHT pairing, not just 'some GEX codes'
        correct = sum(1 for a in atac_bcs if expected[a] in translated)
        print(f"VERDICT: PASS")
        print(f"  --barcode-translate rewrote ATAC barcodes to GEX barcodes.")
        print(f"  {correct}/{len(expected)} expected pairings present in the output.")
        print("  Format confirmed: 2-column TSV written as GEX<TAB>ATAC, i.e.")
        print("  target first, source second (chromap keys on column 2).")
        sys.exit(0)

    if tr_atac > tr_gex:
        print("VERDICT: FAIL - output still carries ATAC barcodes.")
        print("  chromap accepted the file but did not translate. Either the")
        print("  column order is reversed, or the format is not a 2-column TSV.")
        print("  ACTION: try swapping the columns in resolve_arc_translation.py")
        print("  and re-run this test before trusting any multiome output.")
        sys.exit(1)

    print("VERDICT: INCONCLUSIVE")
    print(f"  atac_matches={tr_atac} gex_matches={tr_gex} of {len(expected)} expected")
    sys.exit(1)


if __name__ == "__main__":
    main()

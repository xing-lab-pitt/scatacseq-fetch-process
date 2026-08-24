#!/usr/bin/env python3
"""Decide empirically whether the barcode whitelist needs reverse-complementing.

THE PROBLEM: whether a 10x barcode whitelist matches the reads as written or
reverse-complemented depends on the SEQUENCER's index-read chemistry (NovaSeq
v1.5 / NextSeq / MiSeq differ), not on the assay. It therefore cannot be
hardcoded per protocol, and asking a user to guess produces a silent failure --
chromap runs happily to completion and emits an almost-empty fragment file.

THE FIX (ported from ENCODE_scatac's barcode_revcomp_detect.py): sample reads
from the barcode FASTQ, count exact whitelist matches in BOTH orientations, and
take the winner. If NEITHER orientation clears --min-match-rate, abort here --
6 CPU-seconds instead of 6 CPU-hours to learn the same thing.

The `offset` handling is also from ENCODE: plain 10x reads carry the 16 bp cell
barcode at position 0, but 10x Multiome ATAC carries it at offset 8. Multiome is
out of scope for this pass, but the parameter is kept so M5 is an addition
rather than a rewrite.
"""
import argparse
import gzip
import sys
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BARCODE_LEN = 16
# Plain 10x scATAC: barcode at the start of the read. 10x Multiome ATAC: offset 8.
MODALITY_OFFSET = {"10x": 0, "multiome": 8}


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_whitelist(path):
    with _open(path) as fh:
        return {line.strip().split()[0] for line in fh if line.strip()}


def sample_barcodes(fastq, n_reads, offset):
    """Yield up to n_reads barcode substrings from the FASTQ."""
    out = []
    with _open(fastq) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seq = line.strip()
                out.append(seq[offset:offset + BARCODE_LEN])
                if len(out) >= n_reads:
                    break
    return out


def detect(fastq, whitelist_path, modality, n_reads):
    offset = MODALITY_OFFSET.get(modality, 0)
    whitelist = read_whitelist(whitelist_path)
    if not whitelist:
        raise ValueError(f"whitelist {whitelist_path} is empty")

    barcodes = sample_barcodes(fastq, n_reads, offset)
    if not barcodes:
        raise ValueError(f"{fastq} contains no reads")

    fwd = sum(1 for b in barcodes if b in whitelist)
    rev = sum(1 for b in barcodes if revcomp(b) in whitelist)
    return fwd / len(barcodes), rev / len(barcodes), len(barcodes), offset


def probe_all(fastq, whitelists, n_reads):
    """Match rate for every (modality, orientation) combination.

    WHY THIS EXISTS: modality cannot be read reliably from SRA metadata. Tested
    on GSE219015 (a real multiome study): library names are bare GSM ids with no
    'multiome'/'arc' string, and every library has its OWN BioSample, so nothing
    links the ATAC library to its GEX partner. A name or metadata heuristic
    silently mislabels the whole study.

    But the two chemistries ARE distinguishable from the reads themselves: they
    use different whitelists AND a different barcode offset (0 vs 8). So probe
    all four combinations and let the data decide, exactly as orientation is
    already decided. Returns {(modality, orientation): rate}.
    """
    out = {}
    for modality, wl_path in whitelists.items():
        if not wl_path or not Path(wl_path).exists():
            continue
        offset = MODALITY_OFFSET.get(modality, 0)
        whitelist = read_whitelist(wl_path)
        barcodes = sample_barcodes(fastq, n_reads, offset)
        if not barcodes or not whitelist:
            continue
        out[(modality, "forward")] = sum(
            1 for b in barcodes if b in whitelist) / len(barcodes)
        out[(modality, "revcomp")] = sum(
            1 for b in barcodes if revcomp(b) in whitelist) / len(barcodes)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fastq", required=True, help="barcode read FASTQ (R2)")
    ap.add_argument("--whitelist", required=True, help="10x barcode whitelist")
    ap.add_argument("--modality", default="10x",
                    choices=sorted(MODALITY_OFFSET) + ["unknown"],
                    help="10x or multiome to assert and have checked against the "
                         "reads; 'unknown' to let the reads decide. 'unknown' is "
                         "what prepare_runs.py writes when the metadata does not "
                         "say, instead of guessing.")
    ap.add_argument("--whitelist-alt", default="",
                    help="the OTHER modality's whitelist. When given, all four "
                         "(modality x orientation) combinations are probed and "
                         "reported, so a mislabelled modality is caught by the "
                         "data rather than trusted from the sample sheet.")
    ap.add_argument("--orientation", default="auto",
                    choices=["auto", "forward", "revcomp"],
                    help="auto = detect; forward/revcomp = force (still reported)")
    ap.add_argument("--min-match-rate", type=float, default=0.25)
    ap.add_argument("--n-reads", type=int, default=100000)
    ap.add_argument("--out-whitelist", required=True,
                    help="whitelist written in the orientation the reads use")
    ap.add_argument("--report", required=True)
    ap.add_argument("--translation-in", default="",
                    help="multiome only: ARC translation table (GEX<TAB>ATAC)")
    ap.add_argument("--translation-out", default="",
                    help="multiome only: translation table rewritten into the "
                         "SAME orientation as the whitelist above")
    args = ap.parse_args()

    # modality=unknown: the metadata did not say, so measure it. Probe both
    # chemistries and adopt whichever the reads match, then carry on exactly as
    # if the sheet had named it. There is nothing to disagree with, so the
    # mismatch check below is skipped.
    if args.modality == "unknown":
        if not args.whitelist_alt:
            print("modality=unknown needs --whitelist-alt so both chemistries "
                  "can be probed.", file=sys.stderr)
            sys.exit(1)
        rates = probe_all(args.fastq,
                          {"10x": args.whitelist, "multiome": args.whitelist_alt},
                          args.n_reads)
        print("modality=unknown; measured match rates:")
        for (m, o), rate in sorted(rates.items(), key=lambda kv: -kv[1]):
            print(f"  {m:9} {o:8} {rate:.4f}")
        (best_m, _), best_rate = max(rates.items(), key=lambda kv: kv[1])
        if best_rate < args.min_match_rate:
            print(f"\nNo chemistry matches: best is {best_m} at {best_rate:.4f}, "
                  f"below --min-match-rate {args.min_match_rate}.\n"
                  "Neither whitelist fits these reads. Check that the barcode "
                  "read is really R2 and that the whitelists are the right ones.",
                  file=sys.stderr)
            sys.exit(1)
        # Swap BOTH paths together. Setting only --whitelist would leave
        # --whitelist-alt pointing at the same file, and the cross-modality
        # probe below would then compare a whitelist against itself.
        if best_m != "10x":
            args.whitelist, args.whitelist_alt = args.whitelist_alt, args.whitelist
        args.modality = best_m
        print(f"measured modality={best_m} (match {best_rate:.4f})\n")

    try:
        fwd_rate, rev_rate, n, offset = detect(
            args.fastq, args.whitelist, args.modality, args.n_reads)
    except ValueError as e:
        print(f"Barcode orientation detection failed: {e}", file=sys.stderr)
        sys.exit(1)

    if args.orientation == "auto":
        chosen = "forward" if fwd_rate >= rev_rate else "revcomp"
    else:
        chosen = args.orientation

    best = max(fwd_rate, rev_rate)
    lines = [
        f"fastq\t{args.fastq}",
        f"whitelist\t{args.whitelist}",
        f"modality\t{args.modality}",
        f"barcode_offset\t{offset}",
        f"reads_sampled\t{n}",
        f"forward_match_rate\t{fwd_rate:.4f}",
        f"revcomp_match_rate\t{rev_rate:.4f}",
        f"orientation_mode\t{args.orientation}",
        f"orientation_chosen\t{chosen}",
    ]
    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    Path(args.report).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    # --- cross-modality probe: is the sheet's modality even right? ----------
    if args.whitelist_alt:
        other = "multiome" if args.modality == "10x" else "10x"
        rates = probe_all(args.fastq,
                          {args.modality: args.whitelist, other: args.whitelist_alt},
                          args.n_reads)
        print("\nall (modality, orientation) match rates:")
        for (m, o), rate in sorted(rates.items(), key=lambda kv: -kv[1]):
            mark = "  <- sheet says this modality" if m == args.modality else ""
            print(f"  {m:9} {o:8} {rate:.4f}{mark}")

        if rates:
            (best_m, best_o), best_rate = max(rates.items(), key=lambda kv: kv[1])
            if best_m != args.modality and best_rate > 2 * max(best, 1e-9) \
                    and best_rate >= args.min_match_rate:
                print(
                    f"\nMODALITY MISMATCH: samples.tsv says modality="
                    f"{args.modality}, but the reads match {best_m} much better "
                    f"({best_rate:.4f} vs {best:.4f}).\n"
                    f"Set modality={best_m} for this sample in samples.tsv and "
                    "rerun.\n"
                    "Why this is not auto-corrected: modality also selects the "
                    "barcode-translation behaviour and the downstream pairing, so "
                    "it is a sheet-level decision that should be explicit and "
                    "reviewable, not silently rewritten mid-run.\n"
                    "NOTE: modality is NOT reliably inferable from SRA metadata "
                    "(verified on GSE219015), which is why this check exists.",
                    file=sys.stderr)
                sys.exit(1)

    # Hard stop: neither orientation matches. Downstream this would look like a
    # successful run producing near-zero fragments, so fail here with the numbers.
    if best < args.min_match_rate:
        print(
            f"\nNeither orientation matches the whitelist: forward={fwd_rate:.4f}, "
            f"revcomp={rev_rate:.4f}, both below --min-match-rate="
            f"{args.min_match_rate}.\n"
            "Likely causes, in order of frequency:\n"
            "  1. Wrong whitelist for the chemistry (10x scATAC needs "
            "737K-cratac-v1.txt; 737K-arc-v1.txt is Multiome only).\n"
            "  2. The barcode read was mis-assigned -- check that R2 really is "
            "the 16 bp read (see select_reads.py / prepare_runs.py).\n"
            "  3. The data is not 10x droplet scATAC at all.",
            file=sys.stderr)
        sys.exit(1)

    whitelist = read_whitelist(args.whitelist)
    Path(args.out_whitelist).parent.mkdir(parents=True, exist_ok=True)
    oriented = [bc if chosen == "forward" else revcomp(bc) for bc in sorted(whitelist)]
    with open(args.out_whitelist, "w") as fh:
        fh.write("\n".join(oriented) + "\n")
    print(f"\nWrote {chosen} whitelist ({len(whitelist):,} barcodes) -> "
          f"{args.out_whitelist}")

    # --- multiome: the translation table must match the whitelist orientation --
    #
    # THE BUG THIS FIXES: chromap corrects each read's barcode against the
    # ORIENTED whitelist, then looks that corrected barcode up in the translation
    # table. If the reads are reverse-complemented, the corrected barcodes are
    # RC'd -- but the table built from the ARC whitelists keys on FORWARD ATAC
    # barcodes. chromap then aborts with:
    #     "Barcode does not exist in the translation table."
    # which names neither the barcode nor the cause.
    #
    # So when the orientation is revcomp we rewrite the table's SOURCE column
    # (column 2 -- chromap keys on column 2 and emits column 1) into the same
    # orientation. The target column is left alone: the GEX barcode we want to
    # emit is defined by the GEX whitelist, not by this library's read chemistry.
    if args.translation_out and not args.translation_in:
        # Plain 10x: no translation happens, but the file is still declared as a
        # rule output so the DAG does not need a conditional output (Snakemake
        # handles per-wildcard optional outputs poorly). An empty file makes
        # "no translation" explicit rather than absent.
        Path(args.translation_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.translation_out).write_text("")
        print(f"modality={args.modality}: no barcode translation needed; "
              f"wrote empty {args.translation_out}")

    if args.translation_in:
        if not args.translation_out:
            print("--translation-in given without --translation-out", file=sys.stderr)
            sys.exit(1)
        oriented_set = set(oriented)
        Path(args.translation_out).parent.mkdir(parents=True, exist_ok=True)
        n, matched = 0, 0
        with open(args.translation_in) as fi, open(args.translation_out, "w") as fo:
            for line in fi:
                line = line.rstrip("\n")
                if not line:
                    continue
                to_bc, _, from_bc = line.partition("\t")
                if not from_bc:
                    to_bc, _, from_bc = line.partition(",")
                if chosen == "revcomp":
                    from_bc = revcomp(from_bc)
                fo.write(f"{to_bc}\t{from_bc}\n")
                n += 1
                if from_bc in oriented_set:
                    matched += 1

        cover = matched / len(oriented_set) if oriented_set else 0.0
        print(f"Wrote {chosen}-oriented translation table ({n:,} pairs) -> "
              f"{args.translation_out}")
        print(f"  covers {matched:,}/{len(oriented_set):,} whitelist barcodes "
              f"({cover:.1%})")

        # Fail HERE, with a diagnosis, rather than letting chromap abort later
        # with a message that names neither the barcode nor the cause.
        if cover < 0.99:
            print(
                f"\nTranslation table covers only {cover:.1%} of the oriented "
                "whitelist. chromap would abort mid-alignment with 'Barcode does "
                "not exist in the translation table.'\n"
                "Causes, in order of likelihood:\n"
                "  1. the whitelist and the translation table come from different "
                "ARC versions or different files;\n"
                "  2. the barcode whitelist for this sample is the plain scATAC "
                "list (737K-cratac-v1) rather than the ARC ATAC list;\n"
                "  3. the translation table's column order is wrong (it must be "
                "GEX<TAB>ATAC -- chromap keys on column 2).",
                file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()

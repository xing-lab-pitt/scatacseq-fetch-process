#!/usr/bin/env python3
"""Build the per-sample peak x cell .h5ad from a fragment file and a peak set.

This is the object the rest of the world consumes. It carries:
  X    -- provisional peak x cell counts (snapATAC2)
  var  -- peak coordinates split into chrom / start / end
  obs  -- per-barcode QC: n_fragments, frip, tss_enrichment, nucleosome_signal,
          frac_mito, frac_dup, is_cell
  uns  -- sample-level metrics the QC gate reads, plus 'fragments': the path to
          the fragment file, which is the artifact you never discard (peaks are
          contingent on the population called; fragments are not).

METRIC PROVENANCE -- every number here comes from a source that was verified to
produce it, not from a guessed file format:
  n_fragment, frac_dup, frac_mito  snapATAC2 import_fragments
  frip                             snapATAC2 metrics.frip (peak-set based)
  tss_enrichment                   snapATAC2 metrics.tsse (GTF based)
  nucleosome_signal                computed here in one pass over the fragments
  valid_barcode_frac               the barcode-orientation report's match rate

  mapping_rate                     chromap --summary (see below)

  NOTE: chromap's --summary is deliberately NOT parsed for duplicate/barcode
  rates. Those are computed from the fragments themselves, by the tool that
  produced them. The summary is still written and shipped to MultiQC.

  The one exception is mapping_rate. The fragment file only contains reads that
  mapped, so nothing downstream can tell whether it holds 97% of the library or
  0.7% of it. The summary is the only place that count exists. A run of this
  study produced a complete, QC-passing object from 0.7% of its reads, so this
  metric is now gated.

CELL CALLING is the transparent two-threshold rule from config (scATAC-pro's
argument: the cutoffs cannot be standardised across tissues, so they must at
least be explicit and version-controlled). It is recorded as obs['is_cell'];
no barcode is dropped, so a different call can be made downstream without
re-running anything.
"""
import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pysam
import snapatac2 as snap

MONO_NUCLEOSOME = 147   # bp; the canonical nucleosome footprint


def chrom_sizes_from_fasta(fasta):
    """Chromosome sizes taken from the ALIGNMENT reference itself.

    Deliberately not from snapatac2's bundled genome objects: taking them from
    the same FASTA chromap aligned against is what guarantees the fragment
    contig names and the sizes agree. A 'chr1' vs '1' mismatch here is the
    classic way to get an all-zero matrix with no error message.
    """
    fa = Path(fasta)
    if not (fa.with_suffix(fa.suffix + ".fai").exists()):
        pysam.faidx(str(fa))          # builds the .fai next to the FASTA
    with pysam.FastaFile(str(fa)) as fh:
        return dict(zip(fh.references, fh.lengths))


def nucleosome_signal(fragment_file, barcodes):
    """Per-barcode ratio of mono-nucleosomal to sub-nucleosomal fragments.

    snapATAC2's frag_size_distr is aggregate-only, so this is computed directly:
    one streaming pass over the (already coordinate-sorted, bgzipped) fragment
    file, tallying fragment lengths per barcode. Restricted to `barcodes` so the
    memory cost is bounded by the number of cells, not the whole whitelist.
    """
    wanted = set(barcodes)
    short = defaultdict(int)
    long_ = defaultdict(int)
    with gzip.open(fragment_file, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            bc = parts[3]
            if bc not in wanted:
                continue
            try:
                length = int(parts[2]) - int(parts[1])
            except ValueError:
                continue
            if length >= MONO_NUCLEOSOME:
                long_[bc] += 1
            else:
                short[bc] += 1
    # Ratio, with a 1-count floor on the denominator so a barcode with no
    # sub-nucleosomal fragments yields a large finite number rather than inf.
    return np.array([long_[b] / max(short[b], 1) for b in barcodes], dtype=float)


def sample_metrics(obs, have_tsse=True):
    """Sample-level QC metrics, computed over CALLED CELLS only.

    Pure function of the obs frame so it can be unit-tested without a fragment
    file. That matters: the previous version of this logic was wrong for months
    of nothing, and the gate tests did not catch it because they injected
    pre-made metric values and only exercised the THRESHOLDS, never the
    computation.

    WHY CELLS ONLY: in droplet data most barcodes are empty. Taking a median
    over all barcodes measures the ambient background, and the error GROWS WITH
    SEQUENCING DEPTH, because deeper runs push more empty droplets above the
    import floor. Measured on two real datasets:

                        cells/barcodes  median TSSe (all)  median TSSe (cells)
      atac_pbmc_500        451/594  76%       21.45              22.74
      atac_v1_pbmc_5k    5026/63614   8%        1.36              23.37

    Both are good samples. On the all-barcode median the deeper, better-
    sequenced one scores 17x worse and fails a floor of 4.0. A metric that
    tracks depth rather than quality is worse than no metric.
    """
    cell_mask = obs["is_cell"].to_numpy(dtype=bool)
    n_cells = int(cell_mask.sum())
    out = {"n_cells": n_cells,
           "cell_fraction_of_barcodes": float(n_cells / len(obs)) if len(obs) else 0.0}

    if n_cells:
        cells = obs[cell_mask]
        w = cells["n_fragments"].to_numpy(dtype=float)
        f = cells["frip"].to_numpy(dtype=float)
        # Depth-weighted: a cell with 500 fragments should not count as much as
        # one with 50,000.
        out["frac_reads_in_peaks"] = (
            float(np.average(f, weights=w)) if w.sum() else 0.0)
        out["median_tss_enrichment"] = (
            float(np.nanmedian(cells["tss_enrichment"])) if have_tsse else None)
        out["duplicate_rate"] = (float(np.nanmean(cells["frac_dup"]))
                                 if "frac_dup" in cells else None)
    else:
        out["frac_reads_in_peaks"] = 0.0
        out["median_tss_enrichment"] = None
        out["duplicate_rate"] = None

    # Ambient figures kept for diagnosis only -- useful for spotting a failed
    # cell call, but never what the gate judges.
    out["ambient_median_tss_enrichment"] = (
        float(np.nanmedian(obs.loc[~cell_mask, "tss_enrichment"]))
        if have_tsse and (~cell_mask).any() else None)
    return out


def peak_coords(var_names):
    """'chr1:500000-501000' -> (chrom, start, end) columns."""
    chrom, start, end = [], [], []
    for name in var_names:
        try:
            c, span = str(name).rsplit(":", 1)
            s, e = span.split("-")
            chrom.append(c); start.append(int(s)); end.append(int(e))
        except ValueError:
            chrom.append(str(name)); start.append(-1); end.append(-1)
    return chrom, start, end


def mapping_stats(summary_path):
    """Read-level mapping accounting from chromap's --summary.

    The fragment file only holds reads that mapped, so every metric derived from
    it is conditioned on mapping having worked. This is the only place the
    denominator survives, and without it a run that used 0.7% of its reads is
    indistinguishable from one that used all of them.

    Columns are barcode,total,duplicate,unmapped,lowmapq,... one row per barcode.
    Returns an empty dict rather than raising if the file is absent or the
    schema is not what we expect, since a missing metric is caught by the gate.
    """
    path = Path(summary_path)
    if not path.exists():
        return {}

    total = unmapped = lowmapq = 0
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split(",")
        try:
            i_tot, i_un, i_low = (header.index(c) for c in ("total", "unmapped", "lowmapq"))
        except ValueError:
            return {}
        for line in fh:
            f = line.rstrip("\n").split(",")
            if len(f) <= max(i_tot, i_un, i_low):
                continue
            try:
                total += int(f[i_tot]); unmapped += int(f[i_un]); lowmapq += int(f[i_low])
            except ValueError:
                continue

    if total <= 0:
        return {}
    return {
        "n_read_pairs": total,
        "mapping_rate": (total - unmapped) / total,
        "lowmapq_rate": lowmapq / total,
    }


def read_match_rate(report_path):
    """Barcode match rate for the CHOSEN orientation, from the detector report."""
    if not report_path or not Path(report_path).exists():
        return None
    vals = {}
    for line in Path(report_path).read_text().splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            vals[k] = v
    chosen = vals.get("orientation_chosen")
    key = "forward_match_rate" if chosen == "forward" else "revcomp_match_rate"
    try:
        return float(vals[key])
    except (KeyError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--peaks", required=True, help="MACS3 narrowPeak")
    ap.add_argument("--fasta", required=True, help="alignment reference FASTA")
    ap.add_argument("--gtf", required=True, help="GTF for TSS enrichment")
    ap.add_argument("--chromap-summary", default="")
    ap.add_argument("--orientation-report", default="")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--cell-min-unique-frags", type=float, default=1000)
    ap.add_argument("--cell-min-frip", type=float, default=0.20)
    ap.add_argument("--min-fragments-prefilter", type=int, default=100,
                    help="barcodes below this are not even imported. Kept well "
                         "below the cell-calling threshold so cell calling stays "
                         "an inspectable decision, not an import side effect.")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    work = Path(args.output).parent
    work.mkdir(parents=True, exist_ok=True)
    tmp_import = work / f".{args.sample}.import.h5ad"
    tmp_peaks = work / f".{args.sample}.peakmat.h5ad"
    for stale in (tmp_import, tmp_peaks):
        if stale.exists():
            stale.unlink()

    chrom_sizes = chrom_sizes_from_fasta(args.fasta)
    print(f"reference: {len(chrom_sizes)} contigs from {args.fasta}")

    # chromap writes coordinate-sorted fragments, NOT barcode-sorted.
    data = snap.pp.import_fragments(
        args.fragments, chrom_sizes=chrom_sizes, sorted_by_barcode=False,
        min_num_fragments=args.min_fragments_prefilter, file=str(tmp_import),
    )
    print(f"imported: {data.n_obs:,} barcodes above "
          f"{args.min_fragments_prefilter} fragments")
    if data.n_obs == 0:
        print("No barcodes survived import. The fragment file is empty or its "
              "contig names do not match the reference -- check "
              "qc/barcodes/<sample>.orientation.txt first.", file=sys.stderr)
        sys.exit(1)

    # FRiP: the obs column is named after the dict key.
    snap.metrics.frip(data, {"frip": args.peaks})
    try:
        snap.metrics.tsse(data, args.gtf)
        have_tsse = True
    except Exception as e:                                   # noqa: BLE001
        print(f"WARNING: TSS enrichment failed ({e}). The QC gate will report "
              "tss_enrichment as missing rather than silently passing.",
              file=sys.stderr)
        have_tsse = False

    peak_mat = snap.pp.make_peak_matrix(data, peak_file=args.peaks,
                                        file=str(tmp_peaks))
    adata = peak_mat.to_memory()
    obs_extra = data.obs[:].to_pandas()                       # polars -> pandas
    data.close(); peak_mat.close()

    # --- assemble obs under the schema check_h5ad enforces -------------------
    for col in obs_extra.columns:
        if col not in adata.obs.columns:
            adata.obs[col] = obs_extra[col].to_numpy()

    adata.obs = adata.obs.rename(columns={"n_fragment": "n_fragments",
                                          "tsse": "tss_enrichment"})
    if not have_tsse:
        adata.obs["tss_enrichment"] = np.nan

    adata.obs["nucleosome_signal"] = nucleosome_signal(
        args.fragments, list(adata.obs_names))

    adata.obs["is_cell"] = (
        (adata.obs["n_fragments"] >= args.cell_min_unique_frags)
        & (adata.obs["frip"] >= args.cell_min_frip)
    )

    chrom, start, end = peak_coords(adata.var_names)
    adata.var["chrom"] = chrom
    adata.var["start"] = start
    adata.var["end"] = end

    # --- sample-level metrics the QC gate reads ------------------------------
    #
    # ALL COMPUTED OVER CALLED CELLS, NOT ALL BARCODES. This is not a stylistic
    # choice -- computing them over all barcodes makes them measure the ambient
    # background instead of the sample, and the error grows with sequencing
    # depth. Measured on two real datasets:
    #
    #                     cells/barcodes   median TSSe (all)   median TSSe (cells)
    #   atac_pbmc_500        451/594  76%        21.45               22.74
    #   atac_v1_pbmc_5k    5026/63614   8%         1.36               23.37
    #
    # Both samples are good. The 5k one would have FAILED a TSS floor of 4.0 on
    # the all-barcode median, purely because deeper sequencing yields more
    # ambient barcodes above the import floor. The metric was tracking depth,
    # not quality.
    adata.uns["sample"] = args.sample
    adata.uns["fragments"] = str(Path(args.fragments).resolve())
    adata.uns["peaks_file"] = str(Path(args.peaks).resolve())
    adata.uns.update(sample_metrics(adata.obs, have_tsse=have_tsse))
    adata.uns["valid_barcode_frac"] = read_match_rate(args.orientation_report)
    if args.chromap_summary:
        adata.uns["chromap_summary"] = str(Path(args.chromap_summary).resolve())
        adata.uns.update(mapping_stats(args.chromap_summary))
    n_cells = adata.uns["n_cells"]

    adata.write_h5ad(args.output)
    for stale in (tmp_import, tmp_peaks):
        stale.unlink(missing_ok=True)

    print(f"\nWrote {args.output}")
    print(f"  {adata.n_obs:,} barcodes x {adata.n_vars:,} peaks "
          f"(nnz={adata.X.nnz:,})")
    print(f"  called cells: {n_cells:,} "
          f"(>= {args.cell_min_unique_frags:g} fragments and "
          f">= {args.cell_min_frip:g} FRiP)")
    print(f"  bulk FRiP: {adata.uns['frac_reads_in_peaks']:.4f}")
    print(f"  median TSS enrichment: {adata.uns['median_tss_enrichment']}")


if __name__ == "__main__":
    main()

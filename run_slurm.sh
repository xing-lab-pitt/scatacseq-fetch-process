#!/bin/bash
# =============================================================================
# Launch the scATAC-seq pipeline on SLURM (Snakemake submits each rule as a job).
#
#   sbatch run_slurm.sh                 # run the controller itself as a job
#   ./run_slurm.sh -n                   # or dry-run directly on a login node
#
# Extra args are passed through to snakemake (e.g. -n for a dry run, or a
# specific target path to run one rule).
#
# WHERE THE CHROMAP INDEX IS STORED (custom directory):
#   The index location is a CONFIG setting, not a flag here. In config/config.yaml
#   set reference.chromap_index to the path you want the built index kept at:
#     reference:
#       chromap_index: "/path/to/refs/GRCh38.chromap"   # writable path
#   rule resolve_chromap_index builds it there once and symlinks
#   <workdir>/chromap_index -> it; later runs reuse it. Leave it "" to build into
#   the (ephemeral) workdir instead. To force one build without a full run:
#     sbatch run_slurm.sh <workdir>/chromap_index      # <workdir> = config workdir
#
# CONFIGURE FOR YOUR SITE — ONE STEP:
#   cp activate.local.sh.example activate.local.sh   # then edit the paths
#
# This script sources activate.local.sh automatically. That file is
# .gitignore'd, which is what keeps absolute paths OUT of this committed file
# and the repo portable. The settings it provides:
#   SCATAC_VENV       path to a Python venv to activate (snakemake, executor
#                     plugin, multiqc, snapatac2, macs3, pysam). If unset,
#                     snakemake must already be on PATH.
#   SCATAC_EXTRA_PATH colon-separated dirs to prepend to PATH, for binary tools
#                     not in the venv — chromap, FastQC, sra-tools:
#                       /path/to/chromap:/path/to/FastQC:/path/to/sratoolkit/bin
#   SCATAC_PROFILE    Snakemake profile dir (default: profiles/slurm)
#
# Each may still be overridden per-invocation as a plain environment variable;
# the file assigns with := so an explicit value always wins.
#
# Adjust the #SBATCH lines below (partition, time) to your cluster.
# =============================================================================
#SBATCH -J scatac_ctl
#SBATCH --partition=dept_cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 7-00:00:00
#SBATCH --output=scatac_ctl.%j.log

set -euo pipefail

# Resolve the repo root.
#
# NOT simply $(dirname "$BASH_SOURCE"): under `sbatch`, SLURM copies this script
# to a spool directory, so BASH_SOURCE points there and the profile/ and
# workflow/ dirs are nowhere to be found. Found the hard way on the first real
# submission. Resolution order:
#   1. $SCATAC_REPO if you set it explicitly
#   2. $SLURM_SUBMIT_DIR when running under sbatch (where you ran sbatch from)
#   3. this script's own directory, for ./run_slurm.sh on a login node
if [[ -n "${SCATAC_REPO:-}" ]]; then
    REPO_DIR="${SCATAC_REPO}"
elif [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_SUBMIT_DIR:-}" \
        && -f "${SLURM_SUBMIT_DIR}/workflow/Snakefile" ]]; then
    REPO_DIR="${SLURM_SUBMIT_DIR}"
else
    REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

if [[ ! -f "${REPO_DIR}/workflow/Snakefile" ]]; then
    echo "Cannot locate the pipeline repo: ${REPO_DIR} has no workflow/Snakefile." >&2
    echo "Run sbatch from the repo root, or set SCATAC_REPO=/path/to/repo." >&2
    exit 1
fi
cd "$REPO_DIR"

# --- site configuration ------------------------------------------------------
# Absolute paths differ on every machine, so they must NOT live in this file --
# it is the part that gets committed. Instead they live in activate.local.sh,
# which is .gitignore'd. Copy activate.local.sh.example and edit it.
#
# THE FAILURE THIS PREVENTS: these settings used to come only from the shell you
# happened to run `sbatch` from. A submission from a shell that lacked them died
# two seconds later with "MISSING on PATH: snakemake" -- while an identical
# command from an interactive shell worked. A launcher whose success depends on
# ambient environment is one you cannot trust or hand to anyone else.
#
# Precedence: an explicit environment variable always wins, because
# activate.local.sh assigns with := (assign only if unset). So a one-off
#   SCATAC_VENV=/other/venv sbatch run_slurm.sh ...
# still overrides the file.
if [[ -f "${REPO_DIR}/activate.local.sh" ]]; then
    # shellcheck disable=SC1091
    source "${REPO_DIR}/activate.local.sh"
fi

PROFILE="${SCATAC_PROFILE:-profiles/slurm}"

# Activate a venv if one was provided.
if [[ -n "${SCATAC_VENV:-}" ]]; then
    # shellcheck disable=SC1091
    source "${SCATAC_VENV}/bin/activate"
fi

# Prepend any site-specific tool dirs (chromap, FastQC, sra-tools, ...) so they
# propagate to every SLURM job (SLURM exports the submit environment by default).
if [[ -n "${SCATAC_EXTRA_PATH:-}" ]]; then
    export PATH="${SCATAC_EXTRA_PATH}:$PATH"
fi

# Sanity: fail early if a required binary is missing. sra-tools (prefetch) is
# only needed for source=sra runs, so it is checked but non-fatal.
#
# NOTE: no bedtools / bgzip / tabix here on purpose. Fragment compression and
# indexing go through pysam, and the fragments-in-peaks overlap through pyranges,
# so neither CLI is a dependency.
for t in snakemake chromap macs3 fastqc multiqc; do
    if ! command -v "$t" >/dev/null; then
        echo "MISSING on PATH: $t" >&2
        echo >&2
        echo "Most often this means site configuration was not picked up." >&2
        echo "  repo:            ${REPO_DIR}" >&2
        echo "  activate.local.sh: $([[ -f "${REPO_DIR}/activate.local.sh" ]] && echo present || echo ABSENT)" >&2
        echo "  SCATAC_VENV:     ${SCATAC_VENV:-<unset>}" >&2
        echo "  SCATAC_EXTRA_PATH: ${SCATAC_EXTRA_PATH:-<unset>}" >&2
        echo >&2
        echo "Fix: cp activate.local.sh.example activate.local.sh && edit it." >&2
        exit 1
    fi
done
command -v prefetch >/dev/null || echo "NOTE: prefetch not on PATH (only needed for source=sra runs)" >&2

exec snakemake --profile "$PROFILE" "$@"

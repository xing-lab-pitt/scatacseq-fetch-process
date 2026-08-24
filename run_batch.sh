#!/usr/bin/env bash
# =============================================================================
# run_batch.sh — ONE reconcile -> decide -> (re)launch cycle over the manifest.
#
# Keeps mechanics in a script; keeps judgment with the human. NOT a daemon: it
# performs a single pass and exits. Call it repeatedly, watching progress with
# snakemake_status.sh between calls, until the reconciler reports COMPLETE or
# only human-flagged failures remain.
#
#   ./run_batch.sh [MANIFEST]        MANIFEST defaults to config/manifest.tsv
#   DRY_RUN=1 ./run_batch.sh         reconcile + show what WOULD launch
#   MAX_RETRIES=2 ./run_batch.sh     cap auto-relaunches per study (default 2)
#   LEDGER=results/successful_samples.tsv   append-only success logbook
#
# DECISION PER STUDY, from reconcile.py categories and their `action`:
#
#   RERUN categories -- transient, hand back to Snakemake:
#     missing       not run yet, or an upstream step crashed
#     corrupt       .h5ad unreadable (truncated / killed write). QUARANTINED
#                   first, so Snakemake rebuilds rather than trusting it.
#     no_fragments  uns['fragments'] dangling -- the durable artifact is gone,
#                   so this means re-aligning.
#
#   FLAG categories -- a rerun of identical inputs gives an identical result,
#   so a human decides. NEVER auto-relaunched:
#     qc_fail       failed the deterministic gate
#     read_qc_fail  the raw reads themselves are bad
#     empty_matrix  zero non-zero entries -- a wiring fault (usually a chrom
#                   naming mismatch), not something a retry fixes
#
# WHY THE SPLIT MATTERS: auto-retrying a FLAG category burns hours re-deriving
# the same failure. For ATAC that is not a small waste -- one SRA download can
# take most of a day.
#
# MANIFEST FORMAT (tab-separated, '#' comments ignored):
#   accession <TAB> workdir <TAB> samples_tsv <TAB> configfile
# =============================================================================
set -uo pipefail

REPO_DIR="${SCATAC_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$REPO_DIR" || exit 2

MANIFEST="${1:-config/manifest.tsv}"
DRY_RUN="${DRY_RUN:-0}"
MAX_RETRIES="${MAX_RETRIES:-2}"
LEDGER="${LEDGER:-results/successful_samples.tsv}"
STATE_DIR="${STATE_DIR:-.batch_state}"

if [[ ! -f "$MANIFEST" ]]; then
    echo "manifest not found: $MANIFEST" >&2
    echo "copy config/manifest.example.tsv and edit it." >&2
    exit 2
fi

PY="${SCATAC_VENV:+$SCATAC_VENV/bin/}python"
command -v "$PY" >/dev/null 2>&1 || PY=python

mkdir -p "$STATE_DIR" "$(dirname "$LEDGER")"

echo "=== reconcile ==="
"$PY" workflow/scripts/reconcile.py --manifest "$MANIFEST" \
      --json "$STATE_DIR/reconcile.json" --report "$STATE_DIR/reconcile.tsv"
rc=$?
if [[ $rc -eq 0 ]]; then
    echo
    echo "All studies COMPLETE. Nothing to launch."
    exit 0
elif [[ $rc -ne 1 ]]; then
    echo "reconcile.py failed (exit $rc)" >&2
    exit 2
fi

echo
echo "=== decide ==="
# One launch per study that has any RERUN-category sample and is under its cap.
"$PY" - "$STATE_DIR/reconcile.json" "$MANIFEST" "$STATE_DIR" "$MAX_RETRIES" <<'PYEOF' > "$STATE_DIR/to_launch.txt"
import json, sys, pathlib
res = json.load(open(sys.argv[1]))
manifest, state, cap = sys.argv[2], pathlib.Path(sys.argv[3]), int(sys.argv[4])

cfg = {}
for line in open(manifest):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    if f and f[0].strip().lower() == "accession":
        continue                      # uncommented header row
    if len(f) >= 4:
        cfg[f[0]] = f[3]

for study in res:
    acc = study["accession"]
    rerun = [s for s in study["samples"] if s["action"] == "RERUN"]
    flag = [s for s in study["samples"] if s["action"] == "FLAG"]
    if not rerun:
        if flag:
            print(f"# {acc}: {len(flag)} sample(s) need a human; not relaunching",
                  file=sys.stderr)
        continue
    tries_file = state / f"{acc}.tries"
    tries = int(tries_file.read_text()) if tries_file.exists() else 0
    if tries >= cap:
        print(f"# {acc}: hit MAX_RETRIES={cap}; not relaunching", file=sys.stderr)
        continue
    corrupt = [s["h5ad"] for s in rerun if s["category"] in ("corrupt", "no_fragments")]
    print("\t".join([acc, cfg.get(acc, ""), str(tries + 1), ";".join(corrupt)]))
PYEOF

if [[ ! -s "$STATE_DIR/to_launch.txt" ]]; then
    echo "Nothing to relaunch (only human-flagged failures, or retry caps hit)."
    exit 1
fi

echo
echo "=== launch ==="
while IFS=$'\t' read -r acc configfile tries corrupt; do
    [[ -z "$acc" ]] && continue
    echo "-- $acc (attempt $tries)"

    # Quarantine unusable objects so Snakemake rebuilds them rather than
    # treating "file exists" as done.
    if [[ -n "$corrupt" ]]; then
        IFS=';' read -ra bad <<< "$corrupt"
        for f in "${bad[@]}"; do
            [[ -f "$f" ]] || continue
            if [[ "$DRY_RUN" == "1" ]]; then
                echo "   would quarantine $f"
            else
                mv "$f" "$f.quarantined.$(date +%s)"
                echo "   quarantined $f"
            fi
        done
    fi

    if [[ "$DRY_RUN" == "1" ]]; then
        echo "   would: sbatch run_slurm.sh --configfile $configfile"
        continue
    fi
    if [[ -z "$configfile" ]]; then
        echo "   no configfile in the manifest for $acc; skipping" >&2
        continue
    fi
    jid=$(sbatch --parsable --export=ALL run_slurm.sh --configfile "$configfile" 2>&1)
    if [[ $? -eq 0 ]]; then
        echo "   launched: $jid"
        echo "$tries" > "$STATE_DIR/$acc.tries"
    else
        echo "   sbatch failed: $jid" >&2
    fi
done < "$STATE_DIR/to_launch.txt"

# Append newly-DONE samples to the success logbook (append-only, deduped).
"$PY" - "$STATE_DIR/reconcile.json" "$LEDGER" <<'PYEOF'
import json, sys, pathlib, datetime
res = json.load(open(sys.argv[1])); led = pathlib.Path(sys.argv[2])
seen = set()
if led.exists():
    seen = {l.split("\t")[1] for l in led.read_text().splitlines()[1:] if "\t" in l}
new = [(s["accession"], s["sample"]) for st in res for s in st["samples"]
       if not s["category"] and s["sample"] not in seen]
if new:
    if not led.exists():
        led.write_text("date\tsample\taccession\n")
    with open(led, "a") as fh:
        for acc, samp in new:
            fh.write(f"{datetime.date.today()}\t{samp}\t{acc}\n")
    print(f"ledger: +{len(new)} sample(s) -> {led}")
PYEOF

echo
echo "One cycle done. Re-run this script after the launched jobs finish."
exit 1

#!/usr/bin/env bash
# snakemake_status.sh — DB-free status of a Snakemake + SLURM run.
# Needs only squeue + the log files; does NOT use sacct or seff.
#
# WHY DB-FREE: the SLURM accounting database (slurmdbd) is exactly what goes
# down on this cluster, and when it does, sacct/seff fail while the jobs
# themselves are fine. A status tool that depends on the thing most likely to be
# broken is useless at the moment you need it. This reads only the controller
# log and squeue, which query the controller directly.
#
# Classifies the run as RUNNING / SUCCESS / FAILED / INCOMPLETE from the
# controller log's terminal markers, and on failure prints the failing rule(s)
# and tails their real stderr (the SLURM per-rule logs, not the controller's
# summary -- the controller only says a rule failed, the rule log says why).
#
# Usage:  PIPE=/path/to/repo ./snakemake_status.sh [CONTROLLER_LOG]
#         PIPE defaults to $PWD; the newest scatac_ctl.*.log is used if no log
#         is named.
set -uo pipefail

PIPE="${PIPE:-$PWD}"
LOG="${1:-}"

if [[ -z "$LOG" ]]; then
    mapfile -t _logs < <(ls -t "$PIPE"/scatac_ctl.*.log 2>/dev/null)
    LOG="${_logs[0]:-}"
    # With several runs in one repo, "newest" is often NOT the one you meant --
    # a finished run can be newer than the one still going, so the status would
    # read SUCCESS for the wrong job. Say so rather than answer confidently.
    if [[ ${#_logs[@]} -gt 1 ]]; then
        echo "NOTE: ${#_logs[@]} controller logs in $PIPE; using the newest." >&2
        echo "      Pass one explicitly if you meant a different run:" >&2
        for l in "${_logs[@]:0:4}"; do echo "        $l" >&2; done
    fi
fi
if [[ -z "$LOG" || ! -f "$LOG" ]]; then
    echo "no controller log found in $PIPE (looked for scatac_ctl.*.log)" >&2
    exit 2
fi

echo "log: $LOG"

# --- progress ---------------------------------------------------------------
progress=$(grep -oE '[0-9]+ of [0-9]+ steps \([0-9]+%\) done' "$LOG" | tail -1)
[[ -n "$progress" ]] && echo "progress: $progress"

# --- are its jobs still queued/running? -------------------------------------
ctl_job=$(basename "$LOG" | sed -E 's/scatac_ctl\.([0-9]+)\.log/\1/')
alive=$(squeue -u "$USER" -h -o '%i' 2>/dev/null | grep -c . || true)
ctl_alive=$(squeue -j "$ctl_job" -h -o '%T' 2>/dev/null | grep -c . || true)

# --- RETRY-LOOP DETECTION ---------------------------------------------------
# THE FAILURE THIS EXISTS FOR: Snakemake's SLURM executor submits each rule as a
# job that RE-INVOKES snakemake. A rule that fails and retries does so inside
# that nested session, so "Error in rule" and "Trying to restart" land in the
# RULE's log -- never in the controller log. From the controller's view the job
# is simply still running.
#
# Observed: download_fastq failed and retried 4 times over 15 hours, each cycle
# re-downloading 24 GB and re-extracting ~460 GB, while the controller log
# showed nothing but "Job 4 has been submitted". Directory size kept growing, so
# it even looked like progress.
#
# So: scan the rule logs, not just the controller's.
# SCOPE TO THIS RUN. A first version scanned every rule log from the last 7
# days and flagged a PREVIOUS run's failure as a live loop -- a detector that
# cries wolf on history is one you learn to ignore, which defeats the purpose.
# Only rule logs for job ids this controller actually submitted count.
loops=""
this_run_jobs=$(grep -oE 'submitted with SLURM jobid [0-9]+' "$LOG" 2>/dev/null \
                 | awk '{print $NF}' | sort -u)
if [[ -d "$PIPE/.snakemake/slurm_logs" && -n "$this_run_jobs" ]]; then
    while IFS= read -r rl; do
        [[ -f "$rl" ]] || continue
        # skip logs belonging to any other controller run
        jid=$(basename "$rl" .log)
        grep -qx "$jid" <<< "$this_run_jobs" || continue
        # grep -c prints "0" AND exits 1 on no-match, so "|| echo 0" would
        # append a second 0 and break the arithmetic test below.
        n=$(grep -c 'Trying to restart' "$rl" 2>/dev/null || true)
        n=${n:-0}
        if [[ "$n" -ge 2 ]]; then
            rule=$(basename "$(dirname "$(dirname "$rl")")")
            loops+="  $rule: $n restarts  ($rl)"$'\n'
        fi
    done < <(find "$PIPE/.snakemake/slurm_logs" -name '*.log' 2>/dev/null)
fi

if [[ -n "$loops" ]]; then
    echo
    echo "*** RETRY LOOP DETECTED -- a rule is failing and being retried ***"
    printf '%s' "$loops"
    echo "  A rule retrying >=2 times is almost never transient. Read the rule log"
    echo "  above: the real error is there, not in the controller log."
    echo "  Repeated identical failures burn hours; stop the run and diagnose."
fi

# --- terminal markers -------------------------------------------------------
if grep -q 'steps (100%) done' "$LOG" 2>/dev/null; then
    status="SUCCESS"
elif grep -qE 'Exiting because a job execution failed|WorkflowError|^Error ' "$LOG" 2>/dev/null; then
    status="FAILED"
elif [[ "$ctl_alive" -gt 0 ]]; then
    status="RUNNING"
else
    # Controller gone with no terminal marker: killed, cancelled, or OOM.
    status="INCOMPLETE"
fi
if [[ -n "$loops" && "$status" == "RUNNING" ]]; then
    status="RUNNING-BUT-LOOPING"
fi
echo "status: $status   (controller job $ctl_job, $alive of your jobs queued/running)"

# --- on failure, show the rule logs, not the controller's summary ------------
if [[ "$status" == "FAILED" || "$status" == "INCOMPLETE" ]]; then
    echo
    echo "--- failing rule(s) per the controller ---"
    grep -E 'Error in rule' "$LOG" | sed 's/^/  /' | sort -u | head -5

    echo
    echo "--- real stderr from the rule logs ---"
    # The controller names the per-rule log path; follow it.
    grep -oE '/[^ ]*slurm_logs/[^ )]*\.log' "$LOG" | tail -3 | while read -r rl; do
        [[ -f "$rl" ]] || continue
        echo "  == $rl"
        tail -12 "$rl" | sed 's/^/     /'
    done
fi

case "$status" in
    SUCCESS)              exit 0 ;;
    RUNNING)              exit 3 ;;
    RUNNING-BUT-LOOPING)  exit 5 ;;
    FAILED)     exit 1 ;;
    *)          exit 4 ;;
esac

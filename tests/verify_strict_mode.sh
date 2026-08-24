#!/usr/bin/env bash
# Does qc.mode:strict actually STOP a run?
#
# Every real run so far used mode:warn, but the shipped template defaults to
# strict -- so the mode most users get is the one never exercised end to end.
# qc_gate (reporting) and qc_check (enforcement) are separate rules, and it is
# the wiring between them that has never been tested on real data.
#
# Runs against a scratch workdir seeded with the finished h5ad, so the real
# results are untouched.
set -uo pipefail
R=/net/capricorn/home/xing/lul176/skills_agent/scatacseq-fetch-process
SRC=/net/capricorn/home/xing/data/Luan_data/scatac/GSE219015_mm10
W=${SCRATCH_WD:-/tmp/strict_wd}
V=$(sed -n 's/^: "${SCATAC_VENV:=\(.*\)}"/\1/p' "$R/activate.local.sh")
pass=0; fail=0
ck() { if [ "$2" = "$3" ]; then pass=$((pass+1)); echo "  PASS  $1";
       else fail=$((fail+1)); echo "  FAIL  $1 (got '$2' want '$3')"; fi; }

rm -rf "$W"; mkdir -p "$W/h5ad" "$W/qc/chromap" "$W/qc/peaks" "$W/peaks" "$W/fragments"
cp "$SRC/h5ad/GSM8113084.h5ad"                       "$W/h5ad/"
cp "$SRC/qc/chromap/GSM8113084.chromap_summary.csv"  "$W/qc/chromap/" 2>/dev/null

run_gate () {   # run_gate <min-tss> -> writes report+passlist, echoes qc_check rc
    "$V/bin/python" "$R/workflow/scripts/qc_gate.py" \
        "GSM8113084=$W/h5ad/GSM8113084.h5ad" \
        --report "$W/qc/qc_gate.tsv" --passlist "$W/qc/qc_pass.txt" \
        --min-tss-enrichment "$1" >/dev/null 2>&1
}

echo "=== a sample that PASSES (floor 4.0, actual 15.6) ==="
run_gate 4.0
ck "sample on the passlist" "$(grep -c GSM8113084 "$W/qc/qc_pass.txt")" "1"

echo "=== the same sample FAILING (floor 99) ==="
run_gate 99
ck "passlist is empty"      "$(grep -c . "$W/qc/qc_pass.txt")" "0"
ck "report records the reason" \
   "$(grep -c 'median_tss_enrichment=' "$W/qc/qc_gate.tsv")" "1"
grep -m1 'median_tss_enrichment=' "$W/qc/qc_gate.tsv" | awk -F'\t' '{print "        "$NF}'

echo "=== enforcement: run the REAL qc_check shell body ==="
# Extracted verbatim from rule qc_check rather than reimplemented -- a
# hand-written copy would only prove the copy works. {params}/{input} are
# substituted the way Snakemake substitutes them.
qc_check_body () {   # qc_check_body <strict:True|False> <n_expected>
    local strict="$1" n="$2"
    local passlist="$W/qc/qc_pass.txt" report="$W/qc/qc_gate.tsv" out="$W/qc/qc_ok.txt"
    (
      n_pass=$(grep -c . "$passlist" || true)
      if [ "$strict" = "True" ] && [ "$n_pass" -ne "$n" ]; then
        echo "QC gate: only $n_pass/$n samples passed (mode=strict)." >&2
        echo "See $report. Fix or remove failing samples, then rerun." >&2
        exit 1
      fi
      cp "$passlist" "$out"
    )
}

# Guard against drift: if the rule's body changes, this copy must be updated.
# 5 matching lines: n_pass is assigned, tested, and echoed; plus `exit 1` and
# the final `cp`. If this count moves, rule qc_check changed and the copy above
# needs re-syncing.
body_lines=$(sed -n "/^rule qc_check:/,/^# ---/p" "$R/workflow/Snakefile" | grep -c 'n_pass\|exit 1\|cp {input.passlist}')
ck "rule body unchanged since this test was written" "$body_lines" "5"

rm -f "$W/qc/qc_ok.txt"
qc_check_body True 1 2>"$W/strict.err"; rc_strict=$?
ck "strict + failing sample -> non-zero exit" "$rc_strict" "1"
ck "strict wrote no qc_ok.txt"  "$([ -f "$W/qc/qc_ok.txt" ] && echo yes || echo no)" "no"
ck "strict explains why"        "$(grep -c 'mode=strict' "$W/strict.err")" "1"
sed 's/^/        /' "$W/strict.err" | head -2

qc_check_body False 1 2>/dev/null; rc_warn=$?
ck "warn + failing sample -> exit 0"  "$rc_warn" "0"
ck "warn still produced qc_ok.txt"    "$([ -f "$W/qc/qc_ok.txt" ] && echo yes || echo no)" "yes"

run_gate 4.0
rm -f "$W/qc/qc_ok.txt"
qc_check_body True 1 2>/dev/null; rc_ok=$?
ck "strict + passing sample -> exit 0" "$rc_ok" "0"

rm -rf "$W"
echo
echo "passed: $pass   failed: $fail"
[ "$fail" -eq 0 ]

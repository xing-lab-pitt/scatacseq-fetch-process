#!/usr/bin/env bash
# Do the new guards actually fire on REAL data, not just fixtures?
# A guard that has never fired is indistinguishable from one that is disabled --
# which is exactly what the hg38-blacklist-on-mm10 bug was.
set -uo pipefail
R=/net/capricorn/home/xing/lul176/skills_agent/scatacseq-fetch-process
D=/net/capricorn/home/xing/data/Luan_data/scatac/GSE219015_mm10
B=/net/capricorn/home/xing/lul176/reference/blacklists
FAI=/net/capricorn/home/xing/lul176/reference/mm10/mm10.fa.fai
V=/net/capricorn/home/xing/lul176/mskcc/blood_combined/.venv/bin/python
T=$(mktemp -d)
pass=0; fail=0
ck() { if [ "$2" = "$3" ]; then pass=$((pass+1)); echo "  PASS  $1";
       else fail=$((fail+1)); echo "  FAIL  $1 (got '$2', want '$3')"; fi; }

echo "=== 3a. genome-fai guard: hg38 blacklist against mm10 peaks ==="
$V $R/workflow/scripts/filter_peaks.py \
   --peaks $D/peaks/GSM8113084_peaks.narrowPeak \
   --out $T/x.narrowPeak --report $T/x.tsv \
   --blacklist $B/hg38-blacklist.v2.bed \
   --genome-fai $FAI --primary-chroms-only > $T/a.log 2>&1
ck "aborts on wrong-species blacklist" "$?" "1"
ck "names the offending contigs" \
   "$(grep -c 'absent from the genome' $T/a.log)" "1"
grep -m1 'e.g.' $T/a.log | sed 's/^/        /'

echo "=== 3a-control: correct mm10 blacklist must still pass ==="
$V $R/workflow/scripts/filter_peaks.py \
   --peaks $D/peaks/GSM8113084_peaks.narrowPeak \
   --out $T/ok.narrowPeak --report $T/ok.tsv \
   --blacklist $B/mm10-blacklist.v2.bed \
   --genome-fai $FAI --primary-chroms-only > $T/b.log 2>&1
ck "correct blacklist is accepted" "$?" "0"
ck "drops the same 2016 as the real run" \
   "$(awk -F'\t' '$1=="dropped_blacklisted"{print $2}' $T/ok.tsv)" "2016"

echo "=== 3b. min_mapping_rate gate on the real h5ad (actual rate 0.968) ==="
$V $R/workflow/scripts/qc_gate.py \
   "GSM8113084=$D/h5ad/GSM8113084.h5ad" \
   --report $T/g_hi.tsv --passlist $T/p_hi.txt --min-mapping-rate 0.99 >/dev/null 2>&1
ck "FAILS when threshold above actual" \
   "$(awk -F'\t' 'NR>1{print $NF=="" ? "?" : "x"}' $T/g_hi.tsv | head -1; grep -c . $T/p_hi.txt)" "x
0"

$V $R/workflow/scripts/qc_gate.py \
   "GSM8113084=$D/h5ad/GSM8113084.h5ad" \
   --report $T/g_lo.tsv --passlist $T/p_lo.txt --min-mapping-rate 0.50 >/dev/null 2>&1
ck "PASSES at the configured 0.50" "$(grep -c 'GSM8113084' $T/p_lo.txt)" "1"

echo
echo "  gate verdicts:"
awk -F'\t' 'NR>1{printf "    thr=0.99 -> %s %s\n", $2, $NF}' $T/g_hi.tsv 2>/dev/null | head -2
awk -F'\t' 'NR>1{printf "    thr=0.50 -> %s %s\n", $2, $NF}' $T/g_lo.tsv 2>/dev/null | head -2

rm -rf "$T"
echo
echo "passed: $pass   failed: $fail"

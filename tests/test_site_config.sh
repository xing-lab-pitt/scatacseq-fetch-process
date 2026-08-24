#!/usr/bin/env bash
# Verify run_slurm.sh site-config layering: that the launcher picks up
# activate.local.sh on its own, and that an explicit env var still overrides it.
#
# GUARDS a real failure: a submission from a shell without SCATAC_VENV set died
# 2 s later with "MISSING on PATH: snakemake", while the identical command from
# an interactive shell worked. The launcher depended on ambient environment.
#
# Deliberately machine-independent -- no absolute paths, so it passes on any
# site's activate.local.sh. Skips cleanly where that file has not been created.
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
pass=0; fail=0

if [ ! -f "$R/activate.local.sh" ]; then
    echo "SKIP: no activate.local.sh (cp activate.local.sh.example activate.local.sh)"
    exit 0
fi
ck() { if [ "$2" = "$3" ]; then pass=$((pass+1)); else fail=$((fail+1));
       echo "FAIL: $1"; echo "  expected: $3"; echo "  got:      $2"; fi; }

# 1. syntax
bash -n "$R/run_slurm.sh" && ck "run_slurm.sh parses" ok ok || ck "run_slurm.sh parses" bad ok
bash -n "$R/activate.local.sh" && ck "activate.local.sh parses" ok ok || ck "activate.local.sh parses" bad ok
bash -n "$R/activate.local.sh.example" && ck "example parses" ok ok || ck "example parses" bad ok

# 2. sourcing activate.local.sh yields a REAL venv, starting from a clean env.
#    Asserted as a property (the venv exists and has bin/activate) rather than
#    as a literal path, so this test is not tied to one machine.
got=$(env -u SCATAC_VENV bash -c "source '$R/activate.local.sh'; echo \$SCATAC_VENV")
ck "file sets SCATAC_VENV" "$([ -n "$got" ] && echo set || echo unset)" "set"
ck "SCATAC_VENV is a usable venv" \
   "$([ -x "$got/bin/activate" ] || [ -f "$got/bin/activate" ] && echo yes || echo no)" "yes"

# 3. explicit env var WINS over the file (the := contract)
got=$(SCATAC_VENV=/override/venv bash -c "source '$R/activate.local.sh'; echo \$SCATAC_VENV")
ck "env overrides file" "$got" "/override/venv"

# 4. end-to-end: launcher picks up the venv with NO ambient env at all.
#    This is the exact condition under which job 57147244 died 2s after submit.
#    The regression to guard is "MISSING on PATH" -- NOT the exit status, since
#    snakemake itself will legitimately fail here on the nonexistent configfile.
#    (A stub `snakemake` on PATH does not work as an oracle: activating the venv
#    correctly puts the real one first. Asserting on the stub tested the mock.)
STUB=$(mktemp -d)
for t in chromap macs3 fastqc multiqc; do
    printf '#!/bin/sh\ntrue\n' > "$STUB/$t"
done
chmod +x "$STUB"/*
out=$(env -i HOME="$HOME" PATH="$STUB:/usr/bin:/bin" SCATAC_REPO="$R" \
      bash "$R/run_slurm.sh" --configfile config/__nonexistent__.yaml 2>&1)
ck "no MISSING-on-PATH with empty env" \
   "$(echo "$out" | grep -c 'MISSING on PATH')" "0"
case "$out" in
  *"__nonexistent__.yaml"*) ck "failure came from snakemake, not the launcher" ok ok ;;
  *) ck "failure came from snakemake, not the launcher" "unexpected" "config error" ;;
esac

# 5. with the site file hidden, the error must NAME the fix (not just "MISSING")
mv "$R/activate.local.sh" "$R/activate.local.sh.hidden"
out=$(env -i HOME="$HOME" PATH="/usr/bin:/bin" SCATAC_REPO="$R" \
      bash "$R/run_slurm.sh" 2>&1 | grep -c 'activate.local.sh.example')
mv "$R/activate.local.sh.hidden" "$R/activate.local.sh"
ck "error names the remedy" "$([ "$out" -ge 1 ] && echo yes || echo no)" "yes"

rm -rf "$STUB"
echo "passed: $pass   failed: $fail"
[ "$fail" -eq 0 ]

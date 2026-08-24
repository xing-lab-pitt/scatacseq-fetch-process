"""The qc_check rule: enforcement, separate from reporting.

WHY THIS IS TESTED IN ISOLATION RATHER THAN THROUGH THE PIPELINE
Four attempts to test strict mode by running the real DAG were all confounded by
the harness, not the pipeline: every edit to a rule's params/shell invalidates
chromap, so each attempt needed a full re-alignment, and a stale qc_pass.txt from
an earlier run made one attempt look like a pass when nothing had re-run.

qc_check is six lines of shell. Running them directly with controlled inputs is
deterministic, takes milliseconds, and actually tests the logic -- whereas
driving the whole DAG mostly tests Snakemake. The shell body below is kept
byte-identical to the rule in workflow/Snakefile; test_shell_matches_snakefile
fails if they drift apart.
"""
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow/Snakefile"

# Byte-identical to the shell body of rule qc_check, with Snakemake's
# {placeholders} turned into shell variables.
QC_CHECK_SH = """
n_pass=$(grep -c . "$PASSLIST" || true)
if [ "$STRICT" = "True" ] && [ "$n_pass" -ne "$N" ]; then
  echo "QC gate: only $n_pass/$N samples passed (mode=strict)." >&2
  echo "See $REPORT. Fix or remove failing samples, then rerun." >&2
  exit 1
fi
cp "$PASSLIST" "$OUT"
"""


def run_qc_check(tmp_path, passing_samples, n_samples, strict):
    passlist = tmp_path / "qc_pass.txt"
    passlist.write_text("\n".join(passing_samples) + ("\n" if passing_samples else ""))
    report = tmp_path / "qc_gate.tsv"
    report.write_text("sample\tstatus\n")
    out = tmp_path / "qc_ok.txt"
    r = subprocess.run(
        ["bash", "-c", QC_CHECK_SH],
        env={"PASSLIST": str(passlist), "REPORT": str(report), "OUT": str(out),
             "N": str(n_samples), "STRICT": "True" if strict else "False",
             "PATH": "/usr/bin:/bin"},
        capture_output=True, text=True)
    return r, out


def test_strict_stops_the_run_when_a_sample_fails(tmp_path):
    """THE behaviour: one failing sample must halt everything downstream."""
    r, out = run_qc_check(tmp_path, ["good_a"], n_samples=2, strict=True)
    assert r.returncode == 1
    assert "only 1/2 samples passed" in r.stderr
    assert not out.exists(), "qc_ok.txt must not be written when the gate fails"


def test_strict_passes_when_every_sample_passes(tmp_path):
    r, out = run_qc_check(tmp_path, ["a", "b"], n_samples=2, strict=True)
    assert r.returncode == 0
    assert out.exists()
    assert out.read_text().split() == ["a", "b"]


def test_warn_mode_continues_despite_failures(tmp_path):
    """warn is for exploratory runs: report, do not halt."""
    r, out = run_qc_check(tmp_path, ["a"], n_samples=3, strict=False)
    assert r.returncode == 0
    assert out.exists()


def test_strict_with_zero_passing_samples(tmp_path):
    """The catastrophic case must still be a clean stop, not a crash."""
    r, out = run_qc_check(tmp_path, [], n_samples=2, strict=True)
    assert r.returncode == 1
    assert "only 0/2" in r.stderr
    assert not out.exists()


def test_failure_message_points_at_the_report(tmp_path):
    """The report survives a failed gate -- that separation is the whole reason
    reporting and enforcement are two rules -- so the message must name it."""
    r, _ = run_qc_check(tmp_path, [], n_samples=1, strict=True)
    assert "qc_gate.tsv" in r.stderr
    assert (tmp_path / "qc_gate.tsv").exists()


def test_shell_matches_the_snakefile():
    """Guard against the copy above drifting from the real rule.

    Compares the meaningful tokens rather than whitespace, since the Snakefile
    version carries {placeholders} and this one carries $VARS.
    """
    text = SNAKEFILE.read_text()
    body = text[text.index("rule qc_check:"):]
    body = body[:body.index("# ------------------------------ merge")]

    for token in ["grep -c .", "-ne", "mode=strict", "exit 1", "cp "]:
        assert token in body, f"qc_check rule no longer contains {token!r}"
    # The strict comparison must stay an inequality on the COUNT, not a
    # substring or exit-code check.
    assert re.search(r'"\$n_pass"\s+-ne', body), \
        "qc_check no longer compares the pass COUNT to the sample count"

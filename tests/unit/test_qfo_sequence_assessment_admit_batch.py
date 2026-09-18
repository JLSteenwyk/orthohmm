from pathlib import Path
import subprocess

import pytest


SCRIPT = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_sequence_assessment_admit_batch_20260918.sh"


def test_shell_syntax_and_resource_contract():
    subprocess.run(["bash", "-n", str(SCRIPT)], check=True)
    text = SCRIPT.read_text()
    for expected in ("--cpus-per-task=2", "--mem=64G", "--time=04:00:00", "--no-requeue",
                     "git -C \"$EXECUTOR\" diff --exit-code HEAD -- benchmark_tools",
                     'qfo_sequence_pairs_v1/$VARIANT/results.json',
                     '--job "$JOB" --conversion-job "$CONVERSION_JOB"',
                     '--pairs-sha256 "$SHA"',
                     'qfo_sequence_assessment_admission_${VARIANT}_20260918.json'):
        assert expected in text


@pytest.mark.parametrize("args", [
    [], ["/nonexistent", "a" * 40, "all_hits", "1"],
    ["/nonexistent", "bad", "all_hits", "1", "2"],
    ["/nonexistent", "a" * 40, "unknown", "1", "2"],
    ["/nonexistent", "a" * 40, "all_hits", "x", "2"],
    ["/nonexistent", "a" * 40, "all_hits", "1", "x"],
    ["/nonexistent", "a" * 40, "all_hits", "1", "2", "extra"],
])
def test_malformed_invocation_stops_before_executor_lookup(args):
    done = subprocess.run(["bash", str(SCRIPT), *args], capture_output=True, text=True)
    assert done.returncode != 0
    assert "fatal:" not in done.stderr

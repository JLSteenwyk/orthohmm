"""Prepared launch contracts only; never submit or execute scientific drivers."""

import hashlib
import os
from pathlib import Path
import shlex
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[2]
BATCHES = {
    "pairs": ROOT / "benchmark_tools/results/fault_reviewed_native_qfo_pairs08_20261006.sh",
    "assess": ROOT / "benchmark_tools/results/fault_reviewed_native_qfo_assess08_20261006.sh",
}
DIGEST = "0" * 64


def run_batch(stage, args, **slurm):
    env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
    env.update(slurm)
    return subprocess.run(["bash", str(BATCHES[stage]), *args], env=env,
                          capture_output=True, text=True, timeout=5)


@pytest.mark.parametrize("stage", BATCHES)
def test_bash_syntax(stage):
    result = subprocess.run(["bash", "-n", str(BATCHES[stage])],
                            capture_output=True, text=True, timeout=5)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("args", [[], ["x"], ["f" * 63], ["F" * 64], [DIGEST, "extra"]])
def test_conversion_rejects_bad_arguments_before_driver(args):
    result = run_batch("pairs", args)
    assert result.returncode == 2
    assert "Require the independently verified terminal-review SHA256." in result.stderr


@pytest.mark.parametrize("args", [[], [DIGEST], ["x", "24000"],
    ["F" * 64, "24000"], [DIGEST, "0"], [DIGEST, "-1"], [DIGEST, "024000"],
    [DIGEST, "24000x"], [DIGEST, "24000", "extra"]])
def test_assessment_rejects_bad_arguments_before_driver(args):
    result = run_batch("assess", args)
    assert result.returncode == 2
    assert "Require the verified conversion SHA256 and completed conversion job ID." in result.stderr


@pytest.mark.parametrize("stage,args", [("pairs", [DIGEST]), ("assess", [DIGEST, "24000"])])
def test_not_scheduled_refuses_before_driver(stage, args):
    result = run_batch(stage, args)
    assert result.returncode == 1
    assert not result.stdout and not result.stderr


@pytest.mark.parametrize("stage,args,cpus", [
    ("pairs", [DIGEST], "2"), ("assess", [DIGEST, "24000"], "8")])
@pytest.mark.parametrize("job", ["", "22444", "23017"])
def test_missing_or_reused_job_refuses_before_driver(stage, args, cpus, job):
    result = run_batch(stage, args, SLURM_CPUS_PER_TASK=cpus, SLURM_JOB_ID=job)
    assert result.returncode == 1
    assert not result.stdout and not result.stderr


def test_assessment_cannot_reuse_conversion_job():
    result = run_batch("assess", [DIGEST, "24000"],
                       SLURM_CPUS_PER_TASK="8", SLURM_JOB_ID="24000")
    assert result.returncode == 1
    assert not result.stdout and not result.stderr


@pytest.mark.parametrize("stage,cpus,memory,hours,module,digest", [
    ("pairs", "2", "32G", "06:00:00", "prepare_native_factorial_qfo_pairs",
     "472f26a140cf89aed36f8212a31fef7d8f5061204dac2e0ce80c3c7e6746630e"),
    ("assess", "8", "64G", "04:00:00", "run_native_factorial_qfo_assessment",
     "ba051efa4b2dd9f72a4464f56a891b4742fa3a1c6999f8b95a65d0805c77dafe"),
])
def test_original_driver_and_envelope(stage, cpus, memory, hours, module, digest):
    text = BATCHES[stage].read_text()
    directives = dict(shlex.split(line.removeprefix("#SBATCH "))[0].removeprefix("--").split("=", 1)
                      for line in text.splitlines() if line.startswith("#SBATCH --") and "=" in line)
    assert directives == {
        "job-name": "ohmm_qfo_fault_review_" + ("pairs08" if stage == "pairs" else "assess08"),
        "partition": "gpu", "nodelist": "bizon", "nodes": "1", "ntasks": "1",
        "cpus-per-task": cpus, "mem": memory, "time": hours,
    }
    assert "#SBATCH --no-requeue\n" in text
    assert "#SBATCH --dependency" not in text and "#SBATCH --array" not in text
    assert f'-m benchmark_tools.{module} ' in text.replace("\\\n", " ")
    assert digest in text
    assert hashlib.sha256((ROOT / "benchmark_tools" / (module + ".py")).read_bytes()).hexdigest() == digest
    assert "run_review_gated_native_qfo_pairs" not in text
    assert "run_conversion_gated_native_qfo_assessment" not in text
    assert "-u PYTHONUSERBASE" in text and "-u LD_LIBRARY_PATH" in text
    assert "native_factorial_qfo_pairs_22444_fault_reported_v1" in text
    assert "native_factorial_review_py310_20261004/bin/python" in text


def test_converter_cli_uses_real_review_argument_names_and_separate_destination():
    text = BATCHES["pairs"].read_text()
    assert '--terminal-review "$ROOT/benchmarks/work/native_factorial_terminal_review_22444_fault_reported_v1/review.json"' in text
    assert '--terminal-review-sha256 "$1"' in text
    assert '--output-directory "$ROOT/benchmarks/work/native_factorial_qfo_pairs_22444_fault_reported_v1"' in text
    assert "--review " not in text


def test_assessment_uses_observed_digest_and_conversion_producer():
    text = BATCHES["assess"].read_text()
    assert '--pairs-sha256 "$1" --conversion-job "$2"' in text
    assert '--pairs "$ROOT/benchmarks/work/native_factorial_qfo_pairs_22444_fault_reported_v1/results.json"' in text
    assert "--check-only" not in text and "-resume" not in text

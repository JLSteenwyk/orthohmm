import os
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]
BATCH = ROOT / "benchmark_tools/results/allocated_native_review_23902_20261006.sh"


def test_batch_syntax():
    result = subprocess.run(["bash", "-n", str(BATCH)], capture_output=True, text=True, check=False)
    assert result.returncode == 0 and not result.stdout and not result.stderr


@pytest.mark.parametrize("variables", [{}, {"SLURM_CPUS_PER_TASK": "1", "SLURM_JOB_ID": "99999"},
    {"SLURM_CPUS_PER_TASK": "2", "SLURM_JOB_ID": "0"},
    {"SLURM_CPUS_PER_TASK": "2", "SLURM_JOB_ID": "23902"},
    {"SLURM_CPUS_PER_TASK": "2", "SLURM_JOB_ID": "not-a-job"}])
def test_unscheduled_wrong_cpu_or_same_native_job_refused_before_review(variables):
    environment = dict(os.environ)
    for name in ("SLURM_CPUS_PER_TASK", "SLURM_JOB_ID"):
        environment.pop(name, None)
    environment.update(variables)
    result = subprocess.run(["bash", str(BATCH)], env=environment, capture_output=True, text=True, check=False)
    assert result.returncode == 2
    assert not result.stdout
    assert "Require a distinct scheduled two-CPU review job" in result.stderr


def test_review_only_envelope_request_and_scientific_environment_preserved():
    text = BATCH.read_text()
    for directive in ("--nodes=1", "--ntasks=1", "--cpus-per-task=2", "--mem=128G",
                      "--time=06:00:00", "--no-requeue", "--nodelist=bizon"):
        assert "#SBATCH " + directive in text
    assert "review_allocated_native_factorial_attempt" in text
    assert "request_10_allocated_v1.json" in text
    assert "1355ae3b73a9c8496133bd25a6cbf136d149598ccab680672e81f86482338910" in text
    assert "allocated_native_factorial_terminal_review_23902_v1" in text
    assert "native_factorial_review_py310_20261004/bin/python" in text
    assert "-X faulthandler" in text and "PYTHONUNBUFFERED=1" in text
    assert "n >= 128*2**30" in text and "raw_meminfo=raw" in text
    assert "require(n >= 128*2**30" in text and "assert n" not in text
    for name in ("PYTHONHOME", "PYTHONPATH", "PYTHONUSERBASE", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        assert "-u " + name in text
    assert "run_allocated_native_factorial_cost" not in text
    assert "sbatch" not in text and "scontrol" not in text

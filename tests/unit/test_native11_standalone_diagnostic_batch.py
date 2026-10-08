"""Prospective standalone diagnostic contract; never run the scientific CLI."""

import hashlib
import os
from pathlib import Path
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[2]
BATCH = ROOT / "benchmark_tools/results/native11_standalone_diagnostic_20261008_v1.sh"


def test_bash_syntax():
    result = subprocess.run(["bash", "-n", str(BATCH)], capture_output=True, text=True, timeout=5)
    assert result.returncode == 0 and not result.stdout and not result.stderr


@pytest.mark.parametrize("variables,args", [
    ({}, []),
    ({"SLURM_CPUS_PER_TASK": "1", "SLURM_JOB_ID": "99999"}, []),
    ({"SLURM_CPUS_PER_TASK": "2", "SLURM_JOB_ID": "0"}, []),
    ({"SLURM_CPUS_PER_TASK": "2", "SLURM_JOB_ID": "099999"}, []),
    ({"SLURM_CPUS_PER_TASK": "2", "SLURM_JOB_ID": "23985"}, []),
    ({"SLURM_CPUS_PER_TASK": "2", "SLURM_JOB_ID": "23986"}, []),
    ({"SLURM_CPUS_PER_TASK": "2", "SLURM_JOB_ID": "x"}, []),
    ({"SLURM_CPUS_PER_TASK": "2", "SLURM_JOB_ID": "99999"}, ["unexpected"]),
])
def test_early_guard_refuses_unscheduled_or_reused_identity(variables, args):
    environment = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
    environment.update(variables)
    result = subprocess.run(["bash", str(BATCH), *args], env=environment,
                            capture_output=True, text=True, timeout=5)
    assert result.returncode == 2 and not result.stdout
    assert "Require a distinct scheduled two-CPU diagnostic job and no arguments." in result.stderr


def test_original_scientific_sources_stay_pinned():
    pins = {
        "validate_allocated_native_factorial_outputs.py":
            "146f1333715b5ef52ee17de855fad584d1efc91c7a1154edcd269829b5a871be",
        "validate_native_factorial_outputs.py":
            "3357637503f35238f654edea5c4c12bd293f8d20908c141421a2d40278af7c8d",
        "review_allocated_native_factorial_attempt.py":
            "c8707fdb1855a0d73fba537c001b4a355faff094a38a177a05e179148924c929",
    }
    for name, expected in pins.items():
        assert hashlib.sha256((ROOT / "benchmark_tools" / name).read_bytes()).hexdigest() == expected


def test_scope_environment_and_fresh_namespace():
    text = BATCH.read_text()
    directives = [line for line in text.splitlines() if line.startswith("#SBATCH ")]
    assert directives == [
        "#SBATCH --job-name=ohmm_native11_diagnostic", "#SBATCH --partition=gpu",
        "#SBATCH --nodelist=bizon", "#SBATCH --nodes=1", "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=2", "#SBATCH --mem=32G",
        "#SBATCH --time=06:00:00", "#SBATCH --no-requeue",
    ]
    assert "-m benchmark_tools.validate_allocated_native_factorial_outputs" in text
    assert "request_11_allocated_v1.json" in text
    assert "7bf63b80bd5932b9edbd1b2c5ff3fb77f5557e4f6c64e077045d6a50c8d366a1" in text
    assert 'test ! -e "$DEST"' in text and 'test ! -L "$DEST"' in text
    assert 'mkdir "$DEST"' in text and '--output "$DEST/outputs.json"' in text
    assert "native11_standalone_diagnostic_20261008_v1" in text
    assert "native_factorial_review_py310_20261004/bin/python" in text
    assert "-X faulthandler" in text and "PYTHONUNBUFFERED=1" in text
    assert '/usr/bin/time -v -o "$DEST/time.txt"' in text
    assert "require(n >= 32*2**30" in text and "raw_meminfo=raw" in text
    for name in ("PYTHONHOME", "PYTHONPATH", "PYTHONUSERBASE", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        assert "-u " + name in text
    for forbidden in ("review_allocated_native_factorial_attempt", "gc.disable", "gc.set_threshold",
                      "run_allocated_native_factorial_cost", "sbatch", "scontrol", "--check-only"):
        assert forbidden not in text

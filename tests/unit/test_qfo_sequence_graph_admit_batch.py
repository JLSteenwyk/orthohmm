import os
from pathlib import Path
import subprocess

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_sequence_graph_admit_batch_20260918.sh"


def run(arguments, **env):
    environment = {key: value for key, value in os.environ.items() if not key.startswith("SLURM_")}
    environment.update(env)
    return subprocess.run(["bash", str(SCRIPT), *arguments], env=environment, capture_output=True, text=True)


def arguments():
    return ["/nonexistent/plan.json", "a" * 64, "all_hits", "123", "b" * 64]


def test_shell_syntax_and_frozen_resources():
    subprocess.run(["bash", "-n", str(SCRIPT)], check=True)
    text = SCRIPT.read_text()
    for field in ("#SBATCH --cpus-per-task=2", "#SBATCH --mem=192G", "#SBATCH --no-requeue",
                  "COMMIT=f1e21b09c28f270dc3ef2243bdcad86f212b58a0", "--report-sha256", "--plan-sha256"):
        assert field in text


@pytest.mark.parametrize("args", [[], arguments()[:-1], [*arguments(), "extra"]])
def test_argument_count(args):
    result = run(args)
    assert result.returncode == 2
    assert "Require PLAN" in result.stderr


@pytest.mark.parametrize("index,value", [(0,"relative.json"),(1,"bad"),(4,"A" * 64),(3,"0"),(3,"123_0"),(3,"-1")])
def test_invalid_identity(index, value):
    args = arguments()
    args[index] = value
    result = run(args)
    assert result.returncode == 2
    assert "Invalid absolute" in result.stderr


def test_wrong_variant():
    args = arguments()
    args[2] = "best"
    result = run(args)
    assert result.returncode == 2
    assert "Variant must" in result.stderr


@pytest.mark.parametrize("env", [{}, {"SLURM_JOB_ID":"123","SLURM_CPUS_PER_TASK":"1","SLURM_MEM_PER_NODE":"196608"},
                                {"SLURM_JOB_ID":"123","SLURM_CPUS_PER_TASK":"2","SLURM_MEM_PER_NODE":"65536"}])
def test_wrong_allocation(env):
    result = run(arguments(), **env)
    assert result.returncode == 2
    assert "Require scheduled" in result.stderr

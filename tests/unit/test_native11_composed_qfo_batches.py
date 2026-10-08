import os
from pathlib import Path
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[2]
BATCHES = (
    ("conversion", 2, "12:00:00", 5, "prepare_native11_composed_qfo_pairs"),
    ("assessment", 8, "1-02:00:00", 4, "run_native11_composed_qfo_assessment"),
    ("admission", 2, "06:00:00", 5, "admit_native11_composed_qfo_assessment"),
)


@pytest.mark.parametrize("name,cpus,limit,count,module", BATCHES)
def test_resource_envelope_sanitized_runtime_and_exact_entrypoint(name, cpus, limit, count, module):
    path = ROOT / f"benchmark_tools/results/native11_composed_{name}_20261008_v1.sh"
    text = path.read_text()
    for value in (f"#SBATCH --cpus-per-task={cpus}", "#SBATCH --mem=32G",
        f"#SBATCH --time={limit}", "#SBATCH --no-requeue", "#SBATCH --nodelist=bizon",
        "PYTHONNOUSERSITE=1", "PYTHONDONTWRITEBYTECODE=1", "-u LD_PRELOAD",
        "-u LD_LIBRARY_PATH", "-u LD_AUDIT", "OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1",
        f"-m benchmark_tools.{module}", "native_factorial_review_py310_20261004/bin/python"):
        assert value in text
    assert f'"$#" != {count}' in text
    assert "--source-sha256" in text
    assert "--exclusive" not in text
    assert "-resume" not in text
    subprocess.run(["bash", "-n", str(path)], check=True)


@pytest.mark.parametrize("name,cpus,limit,count,module", BATCHES)
def test_unscheduled_invocation_refuses_before_any_worker_runs(name, cpus, limit, count, module):
    path = ROOT / f"benchmark_tools/results/native11_composed_{name}_20261008_v1.sh"
    env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
    arguments = ["/fixture", "a" * 64, "25001", "25002", "b" * 64][:count]
    result = subprocess.run(["bash", str(path), *arguments], env=env, capture_output=True, text=True)
    assert result.returncode == 2
    assert "Require scheduled" in result.stderr
    assert not result.stdout

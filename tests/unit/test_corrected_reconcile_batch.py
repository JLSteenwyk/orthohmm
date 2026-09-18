import hashlib
import json
import os
from pathlib import Path
import subprocess

import pytest


BATCH = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_corrected_factorial_reconcile_batch_20260918.sh"


def test_resources():
    subprocess.run(["bash", "-n", str(BATCH)], check=True)
    for value in ("--array=0-3%1", "--cpus-per-task=32", "--mem=192G",
                  "--time=48:00:00", "--nodelist=bizon", "--no-requeue"):
        assert "#SBATCH " + value in BATCH.read_text()


@pytest.mark.parametrize("index", range(4))
@pytest.mark.parametrize("problem", [None, "missing_input", "wrong_commit", "dirty", "bad_job", "bad_index"])
def test_handoff(tmp_path, index, problem):
    executor, root = tmp_path / "executor", tmp_path / "root"
    source = executor / "benchmark_tools/run_qfo_corrected_factorial_cell.py"
    source.parent.mkdir(parents=True)
    marker = tmp_path / "called.json"
    source.write_text("import json, os, sys\nfrom pathlib import Path\n"
                      f"Path({str(marker)!r}).write_text(json.dumps([sys.argv[1:], "
                      "{k: os.environ[k] for k in ('PYTHONHASHSEED', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')}]))\n")
    env = dict(os.environ, GIT_AUTHOR_NAME="Test", GIT_AUTHOR_EMAIL="test@example.invalid",
               GIT_COMMITTER_NAME="Test", GIT_COMMITTER_EMAIL="test@example.invalid",
               SLURM_ARRAY_TASK_ID=str(index))
    def git(*args):
        return subprocess.check_output(["git", "-C", str(executor), *args], env=env, text=True).strip()
    git("init", "-q")
    git("add", ".")
    git("-c", "commit.gpgsign=false", "commit", "-qm", "fixture")
    commit = git("rev-parse", "HEAD")
    admission = root / "benchmarks/work/qfo_corrected_candidate_admission_20260918.json"
    admission.parent.mkdir(parents=True)
    admission.write_text("{}\n")
    digest = hashlib.sha256(admission.read_bytes()).hexdigest()
    text = BATCH.read_text()
    root_line = next(line for line in text.splitlines() if line.startswith("ROOT="))
    batch = tmp_path / "batch.sh"
    batch.write_text(text.replace(root_line, f"ROOT={root}"))
    job = "123"
    if problem == "missing_input":
        admission.unlink()
    elif problem == "wrong_commit":
        commit = "0" * 40
    elif problem == "dirty":
        source.write_text(source.read_text() + "# changed\n")
    elif problem == "bad_job":
        job = "pending"
    elif problem == "bad_index":
        env["SLURM_ARRAY_TASK_ID"] = "4"
    result = subprocess.run(["bash", str(batch), str(executor), commit, job], env=env,
                            capture_output=True, text=True)
    if problem:
        assert result.returncode != 0
        assert not marker.exists()
    else:
        assert result.returncode == 0, result.stderr
        argv, observed_env = json.loads(marker.read_text())
        assert argv == ["--root", str(root), "--admission", str(admission),
                        "--admission-sha256", digest, "--admission-job", job,
                        "--environment-manifest", str(root / "benchmark_tools/results/publication_variable_native_methods_20260916.json"),
                        "--index", str(index)]
        assert observed_env == {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                                "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}

import hashlib
import json
import os
from pathlib import Path
import subprocess

import pytest


BATCH = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_corrected_replay_batch_20260918.sh"


def test_batch_syntax_and_resources():
    subprocess.run(["bash", "-n", str(BATCH)], check=True)
    text = BATCH.read_text()
    for setting in ("--nodelist=bizon", "--cpus-per-task=32", "--mem=192G",
                    "--time=24:00:00", "--no-requeue"):
        assert f"#SBATCH {setting}" in text


@pytest.mark.parametrize("problem", [None, "missing_args", "commit_format", "commit_mismatch",
                                    "dirty_executor", "missing_plan", "hash_format", "hash_mismatch"])
def test_batch_provenance_and_handoff(tmp_path, problem):
    executor = tmp_path / "executor"
    executor.mkdir()
    source = executor / "benchmark_tools/run_qfo_corrected_replay.py"
    source.parent.mkdir()
    marker = tmp_path / "unexpected_replay"
    source.write_text("import json, os, sys\nfrom pathlib import Path\n"
                      f"Path({str(marker)!r}).write_text(json.dumps({{'argv': sys.argv[1:], "
                      "'environment': {k: os.environ[k] for k in ('PYTHONHASHSEED', "
                      "'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')}}))\n")
    env = dict(os.environ, GIT_AUTHOR_NAME="Test", GIT_AUTHOR_EMAIL="test@example.invalid",
               GIT_COMMITTER_NAME="Test", GIT_COMMITTER_EMAIL="test@example.invalid")
    def git(*args):
        return subprocess.check_output(["git", "-C", str(executor), *args], env=env, text=True).strip()
    git("init", "-q")
    git("add", ".")
    git("-c", "commit.gpgsign=false", "commit", "-qm", "fixture")
    commit = git("rev-parse", "HEAD")
    plan = tmp_path / "plan.json"
    plan.write_text("{}\n")
    digest = hashlib.sha256(plan.read_bytes()).hexdigest()
    if problem == "commit_format":
        commit = "HEAD"
    elif problem == "commit_mismatch":
        commit = "0" * 40
    elif problem == "dirty_executor":
        source.write_text(source.read_text() + "# changed\n")
    elif problem == "missing_plan":
        plan.unlink()
    elif problem == "hash_format":
        digest = "latest"
    elif problem == "hash_mismatch":
        digest = "0" * 64
    args = [] if problem == "missing_args" else [str(executor), commit, str(plan), digest]
    result = subprocess.run(["bash", str(BATCH), *args], capture_output=True, text=True)
    if problem is not None:
        assert result.returncode != 0
        assert not marker.exists()
    else:
        # Only the batch handoff is exercised; this stub does not run inference.
        assert result.returncode == 0, result.stderr
        observed = json.loads(marker.read_text())
        assert observed["argv"][-4:] == ["--plan", str(plan), "--plan-sha256", digest]
        assert observed["environment"] == {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                                           "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}

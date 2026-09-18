import hashlib
import json
import os
from pathlib import Path
import subprocess

import pytest


RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
STAGES = {
    "prepare": ("prepare_qfo_corrected_factorial.py", "benchmarks/work/qfo_corrected_replay_admission_20260918.json",
                "--admission", "--admission-sha256", "--admission-job", "benchmarks/results/qfo_corrected_factorial_v1"),
    "admit": ("admit_qfo_corrected_candidates.py", "benchmarks/results/qfo_corrected_factorial_v1/manifest.json",
              "--manifest", "--manifest-sha256", "--job", "benchmarks/work/qfo_corrected_candidate_admission_20260918.json"),
}


@pytest.mark.parametrize("stage", STAGES)
def test_syntax_resources(stage):
    path = RESULTS / f"qfo_corrected_candidates_{stage}_batch_20260918.sh"
    subprocess.run(["bash", "-n", str(path)], check=True)
    for setting in ("--nodelist=bizon", "--cpus-per-task=2", "--mem=64G", "--time=24:00:00", "--no-requeue"):
        assert "#SBATCH " + setting in path.read_text()


@pytest.mark.parametrize("stage", STAGES)
@pytest.mark.parametrize("problem", [None, "missing_input", "wrong_commit", "bad_commit", "dirty", "bad_job", "missing_args"])
def test_batch_handoff(tmp_path, stage, problem):
    source_name, relative, input_flag, sha_flag, job_flag, destination = STAGES[stage]
    executor, root = tmp_path / "executor", tmp_path / "root"
    source = executor / "benchmark_tools" / source_name
    source.parent.mkdir(parents=True)
    marker = tmp_path / "called.json"
    source.write_text("import json, os, sys\nfrom pathlib import Path\n"
                      f"Path({str(marker)!r}).write_text(json.dumps({{'argv': sys.argv[1:], 'env': "
                      "{k: os.environ[k] for k in ('PYTHONHASHSEED', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')}}))\n")
    env = dict(os.environ, GIT_AUTHOR_NAME="Test", GIT_AUTHOR_EMAIL="test@example.invalid",
               GIT_COMMITTER_NAME="Test", GIT_COMMITTER_EMAIL="test@example.invalid")
    def git(*args):
        return subprocess.check_output(["git", "-C", str(executor), *args], env=env, text=True).strip()
    git("init", "-q")
    git("add", ".")
    git("-c", "commit.gpgsign=false", "commit", "-qm", "fixture")
    commit = git("rev-parse", "HEAD")
    input_path = root / relative
    input_path.parent.mkdir(parents=True)
    input_path.write_text("{}\n")
    digest = hashlib.sha256(input_path.read_bytes()).hexdigest()
    text = (RESULTS / f"qfo_corrected_candidates_{stage}_batch_20260918.sh").read_text()
    root_line = next(line for line in text.splitlines() if line.startswith("ROOT="))
    batch = tmp_path / "batch.sh"
    batch.write_text(text.replace(root_line, f"ROOT={root}"))
    job = "123"
    if problem == "missing_input":
        input_path.unlink()
    elif problem == "wrong_commit":
        commit = "0" * 40
    elif problem == "bad_commit":
        commit = "HEAD"
    elif problem == "dirty":
        source.write_text(source.read_text() + "# changed\n")
    elif problem == "bad_job":
        job = "pending"
    args = [] if problem == "missing_args" else [str(executor), commit, job]
    result = subprocess.run(["bash", str(batch), *args], capture_output=True, text=True)
    if problem:
        assert result.returncode != 0
        assert not marker.exists()
    else:
        assert result.returncode == 0, result.stderr
        observed = json.loads(marker.read_text())
        assert observed["argv"] == ["--root", str(root), input_flag, str(input_path), sha_flag, digest,
                                    job_flag, job, "--output", str(root / destination)]
        assert observed["env"] == {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                                   "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}

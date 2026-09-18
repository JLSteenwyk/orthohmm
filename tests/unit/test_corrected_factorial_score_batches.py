import hashlib
import json
import os
from pathlib import Path
import subprocess

import pytest

from benchmark_tools.bootstrap_qfo_factorial import CELLS


RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
STAGES = {
    "assessment": ("run_qfo_corrected_factorial_assessment.py", "8", "24:00:00"),
    "score_admit": ("admit_qfo_corrected_factorial_assessment.py", "2", "04:00:00"),
}


@pytest.mark.parametrize("stage", STAGES)
def test_resources(stage):
    _, cpus, wall = STAGES[stage]
    path = RESULTS / f"qfo_corrected_factorial_{stage}_batch_20260918.sh"
    subprocess.run(["bash", "-n", str(path)], check=True)
    for value in (f"--cpus-per-task={cpus}", "--mem=64G", f"--time={wall}",
                  "--nodelist=bizon", "--no-requeue"):
        assert "#SBATCH " + value in path.read_text()


@pytest.mark.parametrize("stage", STAGES)
@pytest.mark.parametrize("index", range(8))
@pytest.mark.parametrize("problem", [None, "missing", "wrong_commit", "dirty", "bad_index", "bad_job"])
def test_handoff(tmp_path, stage, index, problem):
    executor, root = tmp_path / "executor", tmp_path / "root"
    source = executor / "benchmark_tools" / STAGES[stage][0]
    source.parent.mkdir(parents=True)
    marker = tmp_path / "called.json"
    source.write_text("import json, sys\nfrom pathlib import Path\n"
                      f"Path({str(marker)!r}).write_text(json.dumps(sys.argv[1:]))\n")
    env = dict(os.environ, GIT_AUTHOR_NAME="Test", GIT_AUTHOR_EMAIL="test@example.invalid",
               GIT_COMMITTER_NAME="Test", GIT_COMMITTER_EMAIL="test@example.invalid")
    def git(*args):
        return subprocess.check_output(["git", "-C", str(executor), *args], env=env, text=True).strip()
    git("init", "-q")
    git("add", ".")
    git("-c", "commit.gpgsign=false", "commit", "-qm", "fixture")
    commit = git("rev-parse", "HEAD")
    manifest = root / "benchmarks/results/qfo_corrected_factorial_pairs_v1" / CELLS[index] / "results.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text('{"index": ' + str(index) + '}\n')
    digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
    if problem == "missing":
        manifest.unlink()
    elif problem == "wrong_commit":
        commit = "0" * 40
    elif problem == "dirty":
        source.write_text(source.read_text() + "# changed\n")
    text = (RESULTS / f"qfo_corrected_factorial_{stage}_batch_20260918.sh").read_text()
    root_line = next(line for line in text.splitlines() if line.startswith("ROOT="))
    batch = tmp_path / "batch.sh"
    batch.write_text(text.replace(root_line, f"ROOT={root}"))
    passed_index = "8" if problem == "bad_index" else str(index)
    job = "pending" if problem == "bad_job" else "123"
    args = [str(executor), commit, passed_index, job]
    if stage == "score_admit":
        args += ["111"]
    result = subprocess.run(["bash", str(batch), *args], env=env, capture_output=True, text=True)
    if problem:
        assert result.returncode != 0
        assert not marker.exists()
    else:
        assert result.returncode == 0, result.stderr
        expected = ["--root", str(root), "--index", str(index)]
        if stage == "assessment":
            expected += ["--pairs-sha256", digest, "--conversion-job", job]
        else:
            expected += ["--job", job, "--conversion-job", "111", "--pairs-sha256", digest,
                         "--output", str(root / f"benchmarks/work/qfo_corrected_factorial_score_admission_{index}_20260918.json")]
        assert json.loads(marker.read_text()) == expected

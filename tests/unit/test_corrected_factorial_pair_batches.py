import hashlib
import json
import os
from pathlib import Path
import subprocess

import pytest


RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


@pytest.mark.parametrize("mode,memory", [("group", "32G"), ("native", "64G")])
def test_resources(mode, memory):
    path = RESULTS / f"qfo_corrected_{mode}_pairs_batch_20260918.sh"
    subprocess.run(["bash", "-n", str(path)], check=True)
    for value in ("--cpus-per-task=2", f"--mem={memory}", "--time=04:00:00",
                  "--nodelist=bizon", "--no-requeue"):
        assert "#SBATCH " + value in path.read_text()


@pytest.mark.parametrize("index", range(8))
@pytest.mark.parametrize("problem", [None, "missing_candidate", "missing_native", "wrong_commit",
                                    "dirty", "wrong_parity", "bad_job"])
def test_handoff(tmp_path, index, problem):
    mode = "native" if index % 2 else "group"
    executor, root = tmp_path / "executor", tmp_path / "root"
    source = executor / f"benchmark_tools/prepare_qfo_corrected_{mode}_pairs.py"
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
    admission = root / "benchmarks/work/qfo_corrected_candidate_admission_20260918.json"
    admission.parent.mkdir(parents=True)
    admission.write_text("{}\n")
    digest = hashlib.sha256(admission.read_bytes()).hexdigest()
    native = admission.parent / f"qfo_corrected_factorial_native_admission_{index // 2}_20260918.json"
    native.write_text('{"cell": ' + str(index) + '}\n')
    native_digest = hashlib.sha256(native.read_bytes()).hexdigest()
    if problem == "missing_candidate":
        admission.unlink()
    elif problem == "missing_native":
        native.unlink()
    elif problem == "wrong_commit":
        commit = "0" * 40
    elif problem == "dirty":
        source.write_text(source.read_text() + "# changed\n")
    text = (RESULTS / f"qfo_corrected_{mode}_pairs_batch_20260918.sh").read_text()
    root_line = next(line for line in text.splitlines() if line.startswith("ROOT="))
    batch = tmp_path / "batch.sh"
    batch.write_text(text.replace(root_line, f"ROOT={root}"))
    passed_index = index ^ 1 if problem == "wrong_parity" else index
    job = "pending" if problem == "bad_job" else "123"
    result = subprocess.run(["bash", str(batch), str(executor), commit, str(passed_index), job],
                            env=env, capture_output=True, text=True)
    fails = problem is not None and not (problem == "missing_native" and mode == "group")
    if fails:
        assert result.returncode != 0
        assert not marker.exists()
    else:
        assert result.returncode == 0, result.stderr
        expected = ["--root", str(root), "--index", str(index)]
        if mode == "group":
            expected += ["--admission", str(admission), "--admission-sha256", digest, "--admission-job", job]
        else:
            expected += ["--candidate-admission", str(admission), "--candidate-sha256", digest,
                         "--candidate-job", job, "--native-admission", str(native), "--native-sha256", native_digest]
        assert json.loads(marker.read_text()) == expected

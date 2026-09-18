import hashlib
import json
import os
from pathlib import Path
import subprocess

import pytest


BATCH = Path(__file__).resolve().parents[2] / "benchmark_tools/results/qfo_corrected_factorial_admit_batch_20260918.sh"


def test_resources():
    subprocess.run(["bash", "-n", str(BATCH)], check=True)
    for value in ("--cpus-per-task=2", "--mem=64G", "--time=04:00:00",
                  "--nodelist=bizon", "--no-requeue"):
        assert "#SBATCH " + value in BATCH.read_text()


@pytest.mark.parametrize("index", range(4))
@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "failed", "running",
                                    "step_only", "bad_raw", "wrong_task", "dirty", "missing_input"])
def test_handoff(tmp_path, index, problem):
    executor, root = tmp_path / "executor", tmp_path / "root"
    source = executor / "benchmark_tools/admit_qfo_corrected_factorial_cell.py"
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
    task, raw = f"123_{index}", str(200 + index)
    row = f"{task}|{raw}|COMPLETED|0:0\n"
    rows = row + f"{task}.batch|{raw}.batch|COMPLETED|0:0\n"
    if problem == "missing":
        rows = ""
    elif problem == "duplicate":
        rows += row
    elif problem == "failed":
        rows = row.replace("COMPLETED|0:0", "FAILED|1:0")
    elif problem == "running":
        rows = row.replace("COMPLETED", "RUNNING")
    elif problem == "step_only":
        rows = f"{task}.batch|{raw}.batch|COMPLETED|0:0\n"
    elif problem == "bad_raw":
        rows = row.replace(f"|{raw}|", f"|{task}|")
    elif problem == "wrong_task":
        rows = row.replace(task, "999_0")
    elif problem == "dirty":
        source.write_text(source.read_text() + "# changed\n")
    elif problem == "missing_input":
        admission.unlink()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    sacct = bin_dir / "sacct"
    sacct.write_text("#!/home/bizon/anaconda3/bin/python\nimport sys\n"
                     f"assert sys.argv[1:] == ['-j', {task!r}, '--parsable2', '--format=JobID%64,JobIDRaw%64,State,ExitCode']\n"
                     f"print({'JobID|JobIDRaw|State|ExitCode' + chr(10) + rows!r}, end='')\n")
    sacct.chmod(0o755)
    env["PATH"] = str(bin_dir) + os.pathsep + env["PATH"]
    text = BATCH.read_text()
    root_line = next(line for line in text.splitlines() if line.startswith("ROOT="))
    batch = tmp_path / "batch.sh"
    batch.write_text(text.replace(root_line, f"ROOT={root}"))
    result = subprocess.run(["bash", str(batch), str(executor), commit, str(index), "123", "111"],
                            env=env, capture_output=True, text=True)
    if problem:
        assert result.returncode != 0
        assert not marker.exists()
    else:
        assert result.returncode == 0, result.stderr
        assert json.loads(marker.read_text()) == ["--root", str(root), "--index", str(index),
            "--job", raw, "--admission", str(admission), "--admission-sha256", digest,
            "--admission-job", "111", "--output",
            str(root / f"benchmarks/work/qfo_corrected_factorial_native_admission_{index}_20260918.json")]

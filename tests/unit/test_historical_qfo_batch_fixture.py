import json
import shlex
import subprocess
import sys

import pytest


HISTORICAL_PYTHON = "/home/bizon/anaconda3/bin/python"


def script(root="/historical/root"):
    return ("#!/bin/bash\n#SBATCH --cpus-per-task=2\nset -euo pipefail\n"
            f"ROOT={root}\nexec {HISTORICAL_PYTHON} -c "
            "'import json,sys; print(json.dumps(sys.argv[1:]))' \"$ROOT\"\n")


def test_copy_changes_only_declared_bindings(tmp_path, stage_historical_qfo_batch):
    source = tmp_path / "source.sh"
    original = script()
    source.write_text(original)
    root = tmp_path / "root path $literal's"
    destination = stage_historical_qfo_batch(source, tmp_path / "copy.sh", root)
    expected = original.replace("ROOT=/historical/root", "ROOT=" + shlex.quote(str(root)))
    expected = expected.replace(HISTORICAL_PYTHON, shlex.quote(sys.executable))
    assert destination.read_text() == expected
    assert source.read_text() == original
    result = subprocess.run(["bash", str(destination)], check=True, capture_output=True, text=True)
    assert json.loads(result.stdout) == [str(root)]


def test_interpreter_path_is_shell_quoted(tmp_path, monkeypatch, stage_historical_qfo_batch):
    executable = tmp_path / "python path $literal's"
    executable.symlink_to(sys.executable)
    monkeypatch.setattr(sys, "executable", str(executable))
    source = tmp_path / "source.sh"
    source.write_text(script())
    root = tmp_path / "root"
    destination = stage_historical_qfo_batch(source, tmp_path / "copy.sh", root)
    result = subprocess.run(["bash", str(destination)], check=True, capture_output=True, text=True)
    assert json.loads(result.stdout) == [str(root)]


@pytest.mark.parametrize("defect", ["missing_root", "two_roots", "missing_python", "extra_python"])
def test_unknown_markers_fail_before_copy(tmp_path, defect, stage_historical_qfo_batch):
    original = script()
    if defect == "missing_root":
        original = original.replace("ROOT=/historical/root\n", "")
    elif defect == "two_roots":
        original += "ROOT=/second/root\n"
    elif defect == "missing_python":
        original = original.replace(HISTORICAL_PYTHON, "/other/python")
    else:
        original += HISTORICAL_PYTHON + " --version\n"
    source = tmp_path / "source.sh"
    source.write_text(original)
    destination = tmp_path / "copy.sh"
    with pytest.raises(ValueError, match="markers changed"):
        stage_historical_qfo_batch(source, destination, tmp_path / "root")
    assert not destination.exists()
    assert source.read_text() == original


def test_existing_copy_is_not_overwritten(tmp_path, stage_historical_qfo_batch):
    source = tmp_path / "source.sh"
    source.write_text(script())
    destination = tmp_path / "copy.sh"
    destination.write_text("retained failure evidence\n")
    with pytest.raises(FileExistsError):
        stage_historical_qfo_batch(source, destination, tmp_path / "root")
    assert destination.read_text() == "retained failure evidence\n"


def test_two_interpreter_calls_are_relocated(tmp_path, stage_historical_qfo_batch):
    original = script().replace("exec " + HISTORICAL_PYTHON,
                                HISTORICAL_PYTHON + " --version\nexec " + HISTORICAL_PYTHON)
    source = tmp_path / "source.sh"
    source.write_text(original)
    destination = stage_historical_qfo_batch(source, tmp_path / "copy.sh", tmp_path / "root",
                                            python_calls=2)
    assert destination.read_text().count(shlex.quote(sys.executable)) == 2
    assert HISTORICAL_PYTHON not in destination.read_text()
    assert source.read_text() == original

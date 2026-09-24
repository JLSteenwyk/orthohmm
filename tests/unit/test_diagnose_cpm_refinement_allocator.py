import copy
import shutil
import subprocess
import sys

import pytest

from benchmark_tools.diagnose_cpm_refinement_allocator import check, debugger_command, environment, record, run, validate_child


def test_symlink_record_preserves_path_and_checks_target_bytes(tmp_path):
    target = tmp_path / "target"
    target.write_text("original")
    link = tmp_path / "link"
    link.symlink_to(target)
    item = record(link)
    assert item["path"] == str(link)
    check(item)
    check(record(target))
    target.write_text("changed")
    with pytest.raises(ValueError, match="identity changed"):
        check(item)


def test_diagnostic_environment(monkeypatch, tmp_path):
    for key in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH"):
        monkeypatch.setenv(key, "bad")
    monkeypatch.setenv("OMP_NUM_THREADS", "32")
    env = environment(tmp_path)
    assert not set(("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH")) & env.keys()
    assert env["PYTHONMALLOC"] == "debug"
    assert env["PYTHONFAULTHANDLER"] == "1"
    assert env["OMP_NUM_THREADS"] == env["OPENBLAS_NUM_THREADS"] == env["MKL_NUM_THREADS"] == "1"
    assert env["PYTHONPATH"] == str(tmp_path)


@pytest.mark.parametrize("problem", [None, "metadata", "bytes", "record", "duplicate"])
def test_child_check(tmp_path, problem):
    path = tmp_path / "output.txt"
    path.write_text("a a\n" if problem == "duplicate" else "a b\n")
    reference = dict(groups=1, genes=2, modules=["frozen"], accuracy_evaluated=False, output=record(path))
    reference["output"]["path"] = "original"
    child = copy.deepcopy(reference)
    child["output"] = record(path)
    if problem == "metadata":
        child["modules"] = ["changed"]
    elif problem == "bytes":
        reference["output"]["sha256"] = "wrong"
    elif problem == "record":
        child["output"]["path"] = "wrong"
    if problem:
        with pytest.raises(ValueError):
            validate_child(child, reference, path, ["a", "b"])
    else:
        assert validate_child(child, reference, path, ["a", "b"]) == dict(genes=2, groups=1)


def test_unscheduled_run_cannot_create_output(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="allocation"):
        run(tmp_path, "fixture")
    assert not list(tmp_path.iterdir())


def test_debugger_is_single_run_and_preserves_scientific_command():
    child = ["python", "-B", "frozen.py", "--mode", "repeat-refinement"]
    command = debugger_command("/usr/bin/gdb", child)
    assert command[command.index("--args") + 1:] == child
    assert command.count("run") == 1
    assert "continue" not in command
    assert "set auto-load off" in command
    assert "set debuginfod enabled off" in command
    assert "set disable-randomization off" in command
    assert "thread apply all bt" in command


@pytest.mark.parametrize("signal_stop", [False, True])
def test_installed_debugger_preserves_failure_and_collects_stack(signal_stop):
    binary = shutil.which("gdb")
    if binary is None:
        pytest.skip("GDB is not installed")
    source = ("import os,signal,resource; resource.setrlimit(resource.RLIMIT_CORE,(0,0)); "
              "os.kill(os.getpid(),signal.SIGSEGV)") if signal_stop else "raise SystemExit(7)"
    child = subprocess.run(debugger_command(binary, [sys.executable, "-c", source]),
                           text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30)
    if signal_stop:
        assert child.returncode != 0
        assert "received signal SIGSEGV" in child.stdout
        assert "#0" in child.stdout
    else:
        assert child.returncode == 7

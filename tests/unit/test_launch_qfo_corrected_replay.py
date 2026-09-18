import json
import sys
from pathlib import Path

import pytest

from benchmark_tools import launch_qfo_corrected_replay as module


def fixture(tmp_path, monkeypatch):
    executor = tmp_path / "executor"
    admitter = tmp_path / "admitter"
    source = executor / "benchmark_tools/prepare_qfo_corrected_replay.py"
    native = admitter / "benchmark_tools/admit_qfo_corrected_high_sensitivity.py"
    admission_path = tmp_path / "benchmarks/work/qfo_corrected_high_sensitivity_admission_20260918.json"
    primary_path = tmp_path / "benchmark_tools/results/qfo_corrected_primary_commands_20260918.json"
    for path in (source, native, admission_path, primary_path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture")
    launcher = tmp_path / "benchmarks/work/publication_qfo_replay_native_v1"
    output = tmp_path / "benchmarks/results/qfo_corrected_checked_replay_v1"
    primary = {"input_directory": str(tmp_path / "input"), "methods": {module.METHOD: {"native_argv": [sys.executable]}}}
    admission = {"source": module.record(native)}
    checkpoint = tmp_path / "checkpoint"
    plan = {"source": module.record(source), "admission": module.record(admission_path),
            "primary_plan": module.record(primary_path), "input_fastas": [],
            "native_command": module.command_for(Path(sys.executable), launcher, output,
                tmp_path / "input", checkpoint, "sha"), "output_root": str(output), "cwd": str(launcher),
            "environment_overrides": {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                                      "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"},
            "expected_stages": module.STAGES}
    monkeypatch.setattr(module, "validate_admission", lambda *args: (checkpoint, "sha"))
    return tmp_path, plan, admission, primary, executor, admitter


def test_frozen_plan(tmp_path, monkeypatch):
    module.validate_plan(*fixture(tmp_path, monkeypatch))


@pytest.mark.parametrize("key,value", [("source", {}), ("admission", {}), ("primary_plan", {}),
    ("native_command", ["unchecked"]), ("output_root", "/other"), ("cwd", "/other"),
    ("environment_overrides", {}), ("expected_stages", [])])
def test_plan_drift(tmp_path, monkeypatch, key, value):
    args = fixture(tmp_path, monkeypatch)
    args[1][key] = value
    with pytest.raises(ValueError):
        module.validate_plan(*args)


@pytest.mark.parametrize("flag,value", [("--cpu", "64"), ("--matrix", "BLOSUM45"),
    ("--cpm-resolution", "0.2"), ("--leiden-seed", "5"), ("--profile-iterations", "2"),
    ("--profile-min-species", "2"), ("--checkpoint-sha256", "changed")])
def test_scientific_settings_drift(tmp_path, monkeypatch, flag, value):
    args = fixture(tmp_path, monkeypatch)
    command = args[1]["native_command"]
    command[command.index(flag) + 1] = value
    with pytest.raises(ValueError, match="scientific"):
        module.validate_plan(*args)


def test_wrong_admission_source(tmp_path, monkeypatch):
    args = fixture(tmp_path, monkeypatch)
    args[2]["source"] = {}
    with pytest.raises(ValueError, match="binding"):
        module.validate_plan(*args)


def test_wrong_python(tmp_path, monkeypatch):
    args = fixture(tmp_path, monkeypatch)
    args[3]["methods"][module.METHOD]["native_argv"] = ["other-python"]
    with pytest.raises(ValueError):
        module.validate_plan(*args)


@pytest.mark.parametrize("symlink", [False, True])
def test_no_implicit_resume(tmp_path, monkeypatch, symlink):
    args = fixture(tmp_path, monkeypatch)
    output = Path(args[1]["output_root"])
    output.parent.mkdir(parents=True, exist_ok=True)
    if symlink:
        output.symlink_to(tmp_path / "missing")
    else:
        output.mkdir()
    with pytest.raises(FileExistsError):
        module.validate_plan(*args)


@pytest.mark.parametrize("state,cpus,memory,node", [("RUNNING", "2", "64G", "bizon"),
    ("COMPLETED", "1", "64G", "bizon"), ("COMPLETED", "2", "32G", "bizon"),
    ("COMPLETED", "2", "64G", "other")])
def test_completed_gate(tmp_path, monkeypatch, state, cpus, memory, node):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        f"JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n123|{state}|0:0|00:01|{node}|{cpus}|{memory}\n")
    with pytest.raises(ValueError):
        module.prepare(tmp_path, 123, 122)
    assert not list(tmp_path.iterdir())


def test_frozen_revision_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "wrong\n")
    with pytest.raises(ValueError, match="executor"):
        module.frozen(tmp_path, module.EXECUTOR)


def launch_fixture(tmp_path, monkeypatch):
    command = [sys.executable, "/frozen/run_qfo_corrected_replay.py", "--plan", "/plan", "--plan-sha256", "sha"]
    report = {"status": "corrected_replay_dispatch_authorized", "command": command,
              "accuracy_evaluated": False, "publication_ready": False}
    monkeypatch.setattr(module, "prepare", lambda *args: report)
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "196608")
    monkeypatch.setattr(module.os, "uname", lambda: type("Host", (), {"nodename": "bizon"})())
    path = tmp_path / "benchmarks/work/qfo_corrected_replay_dispatch_20260918.json"
    path.parent.mkdir(parents=True)
    calls = []
    monkeypatch.setattr(module.os, "execv", lambda exe, argv: calls.append((exe, argv)))
    return report, path, calls


def test_dispatch_uses_frozen_checked_runner(tmp_path, monkeypatch):
    report, path, calls = launch_fixture(tmp_path, monkeypatch)
    module.launch(tmp_path, 122, 121)
    assert calls == [(sys.executable, report["command"])]
    assert json.loads(path.read_text()) == {**report, "job_id": "123"}
    with pytest.raises(FileExistsError):
        module.launch(tmp_path, 122, 121)


def test_check_only_does_not_authorize_or_exec(tmp_path, monkeypatch):
    report, path, calls = launch_fixture(tmp_path, monkeypatch)
    assert module.launch(tmp_path, 122, 121, True) == report
    assert not path.exists() and calls == []


@pytest.mark.parametrize("key,value", [("SLURM_JOB_ID", ""), ("SLURM_CPUS_PER_TASK", "2"), ("SLURM_MEM_PER_NODE", "65536")])
def test_dispatch_allocation(tmp_path, monkeypatch, key, value):
    _, path, calls = launch_fixture(tmp_path, monkeypatch)
    monkeypatch.setenv(key, value)
    with pytest.raises(ValueError, match="scheduled"):
        module.launch(tmp_path, 122, 121)
    assert not path.exists() and calls == []

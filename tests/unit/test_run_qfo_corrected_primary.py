import json
import os
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import run_qfo_corrected_primary as runner


@pytest.fixture
def execution(tmp_path, monkeypatch):
    source = tmp_path / "inputs"
    source.mkdir()
    fasta = source / "species.fasta"
    fasta.write_text(">gene\nACDE\n")
    root = tmp_path / "outputs"
    plan = {"output_root": str(root), "methods": {}}
    for method in runner.METHODS:
        directory = root / method
        plan["methods"][method] = {
            "output": str(directory), "cwd": str(tmp_path),
            "native_argv": ["/bin/true"],
            "copy_inputs_from": str(source), "copy_inputs_to": str(directory / "input"),
        }
    manifest = tmp_path / "plan.json"
    manifest.write_text(json.dumps(plan))
    digest = runner.record(manifest)["sha256"]
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    calls = []

    def verify(plan):
        calls.append(plan)
        return os.environ.copy(), [runner.record(fasta)]

    monkeypatch.setattr(runner, "verify", verify)
    return plan, manifest, digest, calls


@pytest.mark.parametrize("index", [0, 1])
def test_success_requires_separate_admission(execution, index):
    plan, manifest, digest, calls = execution
    runner.run(manifest, digest, index)
    method = runner.METHODS[index]
    status = json.loads((Path(plan["output_root"]) / "execution" / method / "status.json").read_text())
    assert len(calls) == 2
    assert status["status"] == "process_succeeded_pending_native_admission"
    assert status["accuracy_admitted"] is False
    assert status["native_outputs_validated"] is False
    assert status["exit_code"] == 0
    runner.check(status["timing"])
    if index == 1:
        assert len(status["copied_inputs"]) == 1
        runner.check(status["copied_inputs"][0])
    with pytest.raises(FileExistsError):
        runner.run(manifest, digest, index)


def test_failed_native_command_is_recorded(execution):
    plan, manifest, _, calls = execution
    plan["methods"][runner.METHODS[0]]["native_argv"] = ["/bin/false"]
    manifest.write_text(json.dumps(plan))
    with pytest.raises(RuntimeError, match="Native process failed"):
        runner.run(manifest, runner.record(manifest)["sha256"], 0)
    status = json.loads((Path(plan["output_root"]) / "execution" / runner.METHODS[0] / "status.json").read_text())
    assert status["status"] == "failed" and status["exit_code"] == 1
    assert len(calls) == 2


def test_post_run_provenance_failure_is_recorded(execution, monkeypatch):
    plan, manifest, digest, calls = execution
    original = runner.verify

    def verify(plan):
        if calls:
            raise ValueError("Changed input")
        return original(plan)

    monkeypatch.setattr(runner, "verify", verify)
    with pytest.raises(ValueError, match="Changed input"):
        runner.run(manifest, digest, 0)
    status = json.loads((Path(plan["output_root"]) / "execution" / runner.METHODS[0] / "status.json").read_text())
    assert status["status"] == "failed" and status["accuracy_admitted"] is False


@pytest.mark.parametrize("index", [-1, 2, True])
def test_invalid_index(execution, index):
    _, manifest, digest, _ = execution
    with pytest.raises(ValueError):
        runner.run(manifest, digest, index)


def test_wrong_allocation(execution, monkeypatch):
    _, manifest, digest, _ = execution
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    with pytest.raises(ValueError, match="allocation"):
        runner.run(manifest, digest, 0)


@pytest.mark.parametrize("change", [None, "file", "status", "exit_code", "packages"])
def test_child_resolution_comparison(monkeypatch, change):
    runtime = {"source": {"path": "/inspector"}, "python_executable": {"path": "/python"},
               "package_sources": [], "packages": {"orthofinder": "3.1.5"},
               "child_tools": {"diamond": {"status": "resolved", "file": {"sha256": "abc"}, "exit_code": 0}}}
    observed = json.loads(json.dumps(runtime))
    if change == "packages":
        observed["packages"] = {}
    elif change:
        observed["child_tools"]["diamond"][change] = "changed"
    monkeypatch.setattr(runner, "check", lambda item: None)
    monkeypatch.setattr(runner, "record", lambda item: runtime["python_executable"])
    monkeypatch.setattr(runner.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 0, json.dumps(observed)))
    if change:
        with pytest.raises(ValueError, match="Changed OrthoFinder"):
            runner.verify_child_resolution(runtime, {}, "/venv/bin/python", "/frozen")
    else:
        runner.verify_child_resolution(runtime, {}, "/venv/bin/python", "/frozen")

import json
from pathlib import Path

import pytest

from benchmark_tools import run_qfo_corrected_proteinortho as runner


def test_environment_excludes_container_overrides_and_startup_hooks(monkeypatch):
    monkeypatch.setenv("SINGULARITYENV_PATH", "/foreign")
    monkeypatch.setenv("PERL5OPT", "foreign")
    env = runner.native_environment({"observed_environment": {"PATH": "/frozen", "LD_LIBRARY_PATH": "/libs"}})
    assert env["PATH"] == "/frozen" and env["LD_LIBRARY_PATH"] == "/libs"
    assert env["LC_ALL"] == "C"
    assert "PERL5OPT" not in env and "SINGULARITYENV_PATH" not in env


@pytest.fixture
def execution(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    fasta = source / "species.fasta"
    fasta.write_text(">gene\nACDE\n")
    root = tmp_path / "output"
    plan = {"output_root": str(root), "input_directory": str(source), "cwd": str(root / "input"),
            "native_argv": ["/bin/true"], "input_fastas": [runner.record(fasta)],
            "observed_environment": {"PATH": "/usr/bin:/bin"}}
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan))
    snapshot = {"mock": True}
    runtime = {"status": "proteinortho_runtime_frozen_unrun", "plan": runner.record(plan_path),
               "source": runner.record(runner.__file__), "snapshot": snapshot}
    runtime_path = tmp_path / "runtime.json"
    runtime_path.write_text(json.dumps(runtime))
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setattr(runner, "verify_plan", lambda path: plan)
    monkeypatch.setattr(runner, "runtime_snapshot", lambda *args: snapshot)
    return plan, plan_path, runtime_path, runner.record(runtime_path)["sha256"]


def test_success_is_not_accuracy_admission(execution):
    plan, path, runtime, digest = execution
    runner.run(path, runtime, digest)
    report = json.loads((Path(plan["output_root"]) / "execution.json").read_text())
    assert report["status"] == "process_succeeded_pending_native_admission"
    assert report["accuracy_admitted"] is False and report["native_outputs_validated"] is False
    assert len(report["copied_inputs"]) == 1
    runner.check(report["copied_inputs"][0])
    with pytest.raises(FileExistsError):
        runner.run(path, runtime, digest)


def test_failed_native_command_is_preserved(execution):
    plan, path, runtime, digest = execution
    plan["native_argv"] = ["/bin/false"]
    with pytest.raises(RuntimeError, match="Native process failed"):
        runner.run(path, runtime, digest)
    report = json.loads((Path(plan["output_root"]) / "execution.json").read_text())
    assert report["status"] == "failed" and report["exit_code"] == 1


@pytest.mark.parametrize("stage", [1, 2, 3])
def test_runtime_drift_at_each_gate(execution, monkeypatch, stage):
    plan, path, runtime, digest = execution
    calls = []

    def snapshot(*args):
        calls.append(None)
        return {"mock": len(calls) != stage}

    monkeypatch.setattr(runner, "runtime_snapshot", snapshot)
    with pytest.raises(ValueError, match="Runtime drift"):
        runner.run(path, runtime, digest)
    root = Path(plan["output_root"])
    if stage == 1:
        assert not root.exists()
    else:
        assert json.loads((root / "execution.json").read_text())["status"] == "failed"


def test_wrong_allocation(execution, monkeypatch):
    _, path, runtime, digest = execution
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    with pytest.raises(ValueError, match="allocation"):
        runner.run(path, runtime, digest)

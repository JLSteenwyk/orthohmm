import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import run_qfo_corrected_sonic as runner
from benchmark_tools.prepare_qfo_corrected_sonic import DEFAULT_MODE, TOOLS


def test_environment_is_explicit_and_isolates_bytecode(monkeypatch):
    monkeypatch.setenv("CONDA_PREFIX", "/foreign")
    monkeypatch.setenv("LD_PRELOAD", "/foreign.so")
    env = runner.environment({"output_root": "/fresh"}, {"path": "/frozen/bin"})
    assert "CONDA_PREFIX" not in env and "LD_PRELOAD" not in env and "PYTHONPATH" not in env
    assert env["PYTHONPYCACHEPREFIX"] == "/fresh_pycache"
    assert env["PYTHONDONTWRITEBYTECODE"] == "1" and env["PYTHONNOUSERSITE"] == "1"


@pytest.mark.parametrize("changed", [None, "version", "package_inventory", "diamond"])
def test_resolver_identity(changed):
    original = {"status": "read_only_current_resolution_observed", "version": "2.0.9",
                "default_mode": DEFAULT_MODE, "python": {"sha256": "python"}, "package_inventory": {},
                "tools": {k: {"exit_code": 0, "file": {"sha256": k}} for k in TOOLS}}
    observed = copy.deepcopy(original)
    if changed == "diamond":
        observed["tools"][changed]["file"]["sha256"] = "changed"
    elif changed:
        observed[changed] = "changed"
    if changed:
        with pytest.raises(ValueError):
            runner.compare_resolution(original, observed)
    else:
        runner.compare_resolution(original, observed)


@pytest.fixture
def execution(tmp_path, monkeypatch):
    source = tmp_path / "source.fasta"
    source.write_text(">gene\nACDE\n")
    root = tmp_path / "output"
    plan = {"output_root": str(root), "native_argv": ["/bin/true"], "input_fastas": [runner.record(source)]}
    manifest = tmp_path / "plan.json"
    manifest.write_text(json.dumps(plan))
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    calls = []

    def verify(path):
        calls.append(path)
        return plan, {"PATH": "/usr/bin:/bin"}, []

    monkeypatch.setattr(runner, "verify", verify)
    return plan, manifest, calls


def test_success_does_not_admit_native_outputs(execution):
    plan, manifest, calls = execution
    result = runner.run(manifest)
    assert len(calls) == 2
    assert result["status"] == "process_succeeded_pending_native_admission"
    assert result["accuracy_admitted"] is False and result["native_outputs_validated"] is False
    runner.check(result["copied_inputs"][0])
    with pytest.raises(FileExistsError):
        runner.run(manifest)


def test_check_only_leaves_output_absent(execution):
    plan, manifest, _ = execution
    assert runner.run(manifest, True)["status"] == "preflight_passed_no_inference"
    assert not Path(plan["output_root"]).exists()


def test_failure_recorded(execution):
    plan, manifest, _ = execution
    plan["native_argv"] = ["/bin/false"]
    with pytest.raises(RuntimeError, match="Native process failed"):
        runner.run(manifest)
    status = json.loads((Path(plan["output_root"]) / "execution.json").read_text())
    assert status["status"] == "failed" and status["exit_code"] == 1


def test_postflight_change_recorded(execution, monkeypatch):
    plan, manifest, calls = execution
    original = runner.verify

    def verify(path):
        if calls:
            raise ValueError("Changed runtime")
        return original(path)

    monkeypatch.setattr(runner, "verify", verify)
    with pytest.raises(ValueError, match="Changed runtime"):
        runner.run(manifest)
    assert json.loads((Path(plan["output_root"]) / "execution.json").read_text())["status"] == "failed"


def test_wrong_allocation(execution, monkeypatch):
    _, manifest, _ = execution
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    with pytest.raises(ValueError, match="allocation"):
        runner.run(manifest)


@pytest.mark.parametrize("change", [None, "content", "extra", "missing", "manifest"])
def test_corrected_input_validation(tmp_path, change):
    files = []
    for i in range(78):
        path = tmp_path / f"species{i}.fasta"
        path.write_text(f">gene{i}\nACDE\n")
        files.append(runner.record(path))
    stage = tmp_path / "staging_manifest.json"
    stage.write_text(json.dumps({"input_fastas": files}))
    primary = {"input_directory": str(tmp_path), "inputs": [runner.record(stage), *files]}
    if change == "content":
        Path(files[0]["path"]).write_text("changed")
    elif change == "extra":
        (tmp_path / "unexpected.fasta").write_text(">extra\nACDE\n")
    elif change == "missing":
        Path(files[0]["path"]).unlink()
    elif change == "manifest":
        stage.write_text(json.dumps({"input_fastas": files[:-1]}))
    if change:
        with pytest.raises((ValueError, FileNotFoundError)):
            runner.verify_corrected_inputs(primary, files)
    else:
        runner.verify_corrected_inputs(primary, files)

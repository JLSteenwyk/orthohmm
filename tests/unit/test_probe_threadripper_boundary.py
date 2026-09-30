import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import probe_threadripper_boundary as module
from tests.unit.test_replay_threadripper_boundary import archive, evidence


def save(path, value):
    path.write_text(json.dumps(value))


@pytest.fixture
def binding(tmp_path, monkeypatch):
    root = tmp_path / "repository"
    helpers = root / "benchmark_tools"
    helpers.mkdir(parents=True)
    source = helpers / "probe_threadripper_boundary.py"
    source.write_text("# Synthetic source binding, not executed.\n")
    monkeypatch.setattr(module, "__file__", str(source))
    script = root / "submission.sh"
    script.write_text("# Synthetic submission binding, not submitted.\n")
    provenance = root / "controller.json"
    save(provenance, {"synthetic": True})
    policy = root / "policy.md"
    policy.write_text("Synthetic test policy, not execution approval.\n")
    work = root / "benchmarks/work"
    work.mkdir(parents=True)
    path = root / "protocol.json"
    value = dict(schema="threadripper_boundary_fixture_v1", command=["/usr/bin/true"],
        settings=module.SETTINGS, scheduling=module.SCHEDULING, release_guard=None,
        production_execution_authorized=False, **{k: False for k in module.ADMISSIONS},
        cwd=str(root), run_directory=str(work / "fixture"),
        interpreter=module.record(sys.executable), native_executable=module.record("/usr/bin/true"),
        sources=[module.record(source)], submission_script=module.record(script),
        controller_provenance=module.record(provenance), policy=module.record(policy))
    save(path, value)
    return root, path, value


def digest(path):
    return module.record(path)["sha256"]


def relocated_measurement(source, destination):
    destination.mkdir()
    for path in source.iterdir():
        if path.is_file():
            shutil.copy2(path, destination / path.name)
    path = destination / "boundary_report.json"
    measured = json.loads(path.read_text())
    measured["launched"][-1] = str(destination)
    for pin in measured["point_records"]:
        pin["path"] = str(destination / Path(pin["path"]).name)
    save(path, measured)


def environment(root, monkeypatch):
    monkeypatch.chdir(root)
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    for key, value in dict(SLURM_JOB_ID="21816", SLURM_CPUS_PER_TASK="64", SLURM_MEM_PER_NODE="131072").items():
        monkeypatch.setenv(key, value)
    for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        monkeypatch.delenv(key, raising=False)


def test_run_and_independent_raw_audit(binding, archive, monkeypatch):
    root, path, value = binding
    environment(root, monkeypatch)
    directory = Path(value["run_directory"])
    calls = []
    def measure(command, output, job, **kwargs):
        calls.append((command, job, kwargs))
        relocated_measurement(archive, output)
    monkeypatch.setattr(module, "measure", measure)
    result = module.run(directory, path, digest(path))
    assert calls == [(["/usr/bin/true"], 21816, module.SETTINGS)]
    assert result["native_points"] == 2 and result["native_exit_code"] == 0
    assert result["cpu"]["cpu_seconds"] > 0
    assert result["step_memory"]["bytes"] == 200
    assert all(result[k] is False for k in module.ADMISSIONS)
    reread = module.audit(directory, 21816, path, digest(path), directory / "independent.json")
    assert reread == result
    with pytest.raises(FileExistsError):
        module.run(directory, path, digest(path))
    with pytest.raises(FileExistsError):
        module.audit(directory, 21816, path, digest(path), directory / "independent.json")
    with pytest.raises(ValueError):
        module.audit(directory, 21817, path, digest(path), directory / "wrong-job.json")


@pytest.mark.parametrize("fault", ["schema", "command", "settings", "bool_interval", "scheduling",
    "release_guard", "authorized", *module.ADMISSIONS, "cwd", "directory", "sources", "native",
    "pin", "inventory", "digest", "symlink"])
def test_protocol_rejects_changed_bindings(binding, fault):
    root, path, value = binding
    if fault == "schema": value["schema"] = "other"
    elif fault == "command": value["command"] = ["/usr/bin/false"]
    elif fault == "settings": value["settings"] = dict(module.SETTINGS, cpus=20)
    elif fault == "bool_interval": value["settings"] = dict(module.SETTINGS, interval_s=True)
    elif fault == "scheduling": value["scheduling"] = dict(module.SCHEDULING, attempts=2)
    elif fault == "release_guard": value["release_guard"] = "approved"
    elif fault == "authorized": value["production_execution_authorized"] = True
    elif fault in module.ADMISSIONS: value[fault] = True
    elif fault == "cwd": value["cwd"] = str(root.parent)
    elif fault == "directory": value["run_directory"] = str(root / "elsewhere")
    elif fault == "sources": value["sources"] = []
    elif fault == "native": value["native_executable"] = module.record("/usr/bin/false")
    elif fault == "pin": Path(value["submission_script"]["path"]).write_text("changed")
    elif fault == "inventory": (root / "benchmark_tools/extra.py").write_text("changed")
    save(path, value)
    expected = "0" * 64 if fault == "digest" else digest(path)
    if fault == "symlink":
        link = root / "indirect.json"
        link.symlink_to(path)
        path = link
    with pytest.raises(ValueError):
        module.protocol(path, expected)


@pytest.mark.parametrize("fault", ["host", "cpu", "memory", "cwd", "runtime", "preload", "job"])
def test_invalid_execution_environment_does_not_create_output(binding, monkeypatch, fault):
    root, path, value = binding
    environment(root, monkeypatch)
    if fault == "host": monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="dgx"))
    elif fault == "cpu": monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    elif fault == "memory": monkeypatch.setenv("SLURM_MEM_PER_NODE", "65536")
    elif fault == "cwd": monkeypatch.chdir(root.parent)
    elif fault == "runtime":
        value["interpreter"] = module.record("/usr/bin/true")
        save(path, value)
    elif fault == "preload": monkeypatch.setenv("LD_PRELOAD", "untrusted")
    else: monkeypatch.setenv("SLURM_JOB_ID", "0")
    output = Path(value["run_directory"])
    with pytest.raises(ValueError):
        module.run(output, path, digest(path))
    assert not output.exists()


@pytest.mark.parametrize("fault", ["collection", "native_nonzero", "native_timeout", "source_change"])
def test_failed_attempts_are_retained_and_never_retried(binding, archive, monkeypatch, fault):
    root, path, value = binding
    environment(root, monkeypatch)
    output = Path(value["run_directory"])
    calls = []
    def measure(command, directory, job, **kwargs):
        calls.append(job)
        if fault == "collection":
            raise RuntimeError("collection failed")
        relocated_measurement(archive, directory)
        if fault == "source_change":
            Path(module.__file__).write_text("changed")
        else:
            done = json.loads((directory / "done.json").read_text())
            done["exit_code"] = 7 if fault == "native_nonzero" else 124
            done["timed_out"] = fault == "native_timeout"
            save(directory / "done.json", done)
            measured = json.loads((directory / "boundary_report.json").read_text())
            measured.update(native=done, status="command_failed")
            save(directory / "boundary_report.json", measured)
    monkeypatch.setattr(module, "measure", measure)
    with pytest.raises((RuntimeError, ValueError)):
        module.run(output, path, digest(path))
    assert calls == [21816]
    assert (output / "started.json").exists() and (output / "fixture_failure.json").exists()
    assert not (output / "component_audit.json").exists()
    with pytest.raises(ValueError):
        module.audit(output, 21816, path, digest(path), output / "independent.json")

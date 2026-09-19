import json
from pathlib import Path

import pytest

from benchmark_tools import replay_root_context_control as module
from benchmark_tools.root_context_control_design import service_command, unit_name
from tests.unit.test_root_context_control_trial import user_trial, MANAGER, evaluate, save

REMOTE = Path("/deployed/output/trial_03")
RECIPE = Path("/deployed/source")
PYTHON = "/usr/bin/python"


@pytest.fixture
def archive(user_trial, monkeypatch):
    directory, measured, context = user_trial
    work = directory / "workload"
    for kind in ("ready", "done"):
        value = json.loads((work / f"workload_{kind}.json").read_text())
        for row in value["workers"]:
            save(work / f"{kind}_{row['cpu']}.json", row)
    save(work / "workload_go.json", dict(go=True))
    (work / "service.log").touch()
    controller = json.loads((directory / "controller.json").read_text())
    controller.update(unit=unit_name(123, 3), manager=MANAGER,
        command=service_command(PYTHON, RECIPE / "benchmark_tools/full_node_control_workload.py",
                                REMOTE / "workload", unit_name(123, 3), 0))
    save(directory / "controller.json", controller)
    save(directory / "trial.json", evaluate(user_trial))
    calls = []
    monkeypatch.setattr(module, "replay_measurement", lambda *a, **k: calls.append((a, k)) or
        dict(lineage=dict(measured=measured), context=context["context"]))
    return directory, calls


def replay(directory):
    return module.replay(directory, REMOTE, RECIPE, PYTHON, "user-contended", 123, 3, MANAGER)


def test_raw_witness_replay_and_expected_native_command(archive):
    directory, calls = archive
    result = replay(directory)
    assert result["status"] == "root_context_control_replayed"
    assert result["scientific_timings_admitted"] is False
    assert calls[0][0][2] == [PYTHON, "-B", str(RECIPE / "benchmark_tools/full_node_control_workload.py"),
        "--directory", str(REMOTE / "workload"), "--worker", "steady"]
    assert calls[0][1] == dict(expected_timeout_s=60)


@pytest.mark.parametrize("fault", ["raw", "missing", "extra", "go", "command", "manager", "exit_type",
                                  "cleanup", "removal", "summary", "symlink"])
def test_tampered_archive_rejected(archive, fault):
    directory, _ = archive
    work = directory / "workload"
    if fault == "raw":
        path = work / "done_0.json"
        value = json.loads(path.read_text())
        value["self_cpu_s"] = 15
        save(path, value)
    elif fault == "missing":
        (work / "ready_0.json").unlink()
    elif fault == "extra":
        (work / "failed_0.json").touch()
    elif fault == "go":
        save(work / "workload_go.json", dict(go=False))
    elif fault == "summary":
        save(directory / "trial.json", {})
    elif fault == "symlink":
        (work / "link").symlink_to(work / "done_0.json")
    else:
        path = directory / "controller.json"
        value = json.loads(path.read_text())
        if fault == "command":
            value["command"].append("--unexpected")
        elif fault == "manager":
            value["manager"] = "/user.slice/other"
        elif fault == "exit_type":
            value["service_exit_code"] = False
        elif fault == "cleanup":
            value["owned_service_stop"] = dict(returncode=0)
        else:
            value["removal"][0]["monotonic_ns"] = 1
        save(path, value)
    with pytest.raises(ValueError):
        replay(directory)


def test_wrong_frozen_condition_rejected(archive):
    directory, _ = archive
    with pytest.raises(ValueError, match="frozen index"):
        module.replay(directory, REMOTE, RECIPE, PYTHON, "steady", 123, 3, MANAGER)

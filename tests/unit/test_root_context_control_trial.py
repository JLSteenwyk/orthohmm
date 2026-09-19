import json
import threading
from types import SimpleNamespace

import pytest

from benchmark_tools import run_root_context_control_trial as module
from benchmark_tools.root_context_control_design import ORDER, check_service_scope, unit_name
from tests.unit.test_validate_full_node_control import fixture

MANAGER = "/user.slice/user-1000.slice/user@1000.service"


def save(path, value):
    path.write_text(json.dumps(value))


@pytest.fixture
def user_trial(tmp_path):
    work = tmp_path / "workload"
    work.mkdir()
    args = fixture("contended")
    scope = MANAGER + "/app.slice/" + unit_name(123, 3)
    member = "0::" + scope + "\n"
    args["competitor_ready"]["membership"] = member
    args["competitor"].update(membership=member, final_membership=member)
    for name, data in (("workload_ready", args["ready"]), ("workload_done", args["done"]),
                       ("competitor_ready", args["competitor_ready"]), ("competitor_done", args["competitor"])):
        save(work / (name + ".json"), data)
    save(tmp_path / "controller.json", dict(status="completed", service_exit_code=0,
        service_scope=scope, removal=[dict(exists=False, monotonic_ns=24_000_000_000)]))
    points = [dict(native_membership="native", host=[dict(started_monotonic_ns=t,
        raw=dict(cgroup_membership="batch"))], root_context=dict(
            host_before=dict(started_ns=t+1), host_after=dict(finished_ns=t+100)))
        for t in (0, 3_000_000_000, 10_000_000_000, 21_000_000_000, 25_000_000_000)]
    measured = dict(native=args["native"], points=points)
    context = dict(context=dict(observation_window=dict(scope_cpu_usec={"/user.slice": 9_000_000})))
    return tmp_path, measured, context


def evaluate(value):
    directory, measured, context = value
    return module.validate_trial(directory, "user-contended", 123, 3, MANAGER, measured, context)


def test_user_scope_is_explicit_and_common_windows_include_supplement(user_trial):
    result = evaluate(user_trial)
    assert result["positive_control_response"] is True
    assert result["common_intervals"] == [1, 2]
    assert result["scientific_timings_admitted"] is False


def test_low_user_response_remains_false_not_invalid_workload(user_trial):
    user_trial[2]["context"]["observation_window"]["scope_cpu_usec"]["/user.slice"] = 4_999_999
    result = evaluate(user_trial)
    assert result["status"] == "root_context_workload_validated"
    assert result["positive_control_response"] is False


@pytest.mark.parametrize("fault", ["dose", "scope", "removal", "exit", "coordinator", "cleanup", "coverage", "no_common", "native"])
def test_invalid_workload_or_service_evidence_rejected(user_trial, fault):
    directory, measured, context = user_trial
    controller = json.loads((directory / "controller.json").read_text())
    if fault in ("dose", "scope"):
        path = directory / "workload/competitor_done.json"
        data = json.loads(path.read_text())
        if fault == "dose":
            data["self_cpu_s"] = 4
        else:
            data["membership"] = data["final_membership"] = "0::/system.slice/unrelated\n"
        save(path, data)
    elif fault == "removal":
        controller["removal"][-1]["exists"] = True
    elif fault == "exit":
        controller["service_exit_code"] = 1
    elif fault == "coordinator":
        controller["status"] = "failed"
    elif fault == "cleanup":
        controller["cleanup_error"] = "Service stop failed"
    elif fault == "coverage":
        measured["points"][-1]["root_context"]["host_after"]["finished_ns"] = 21_000_000_000
    elif fault == "no_common":
        measured["points"] = [measured["points"][0], measured["points"][-1]]
    else:
        measured["native"]["exit_code"] = 1
    save(directory / "controller.json", controller)
    with pytest.raises(ValueError):
        evaluate(user_trial)


def test_idle_has_no_injected_workload(tmp_path):
    (tmp_path / "workload").mkdir()
    points = [dict(host=[dict(started_monotonic_ns=t)], root_context=dict(
        host_before=dict(started_ns=t+1), host_after=dict(finished_ns=t+2)))
        for t in (1, 2_000_000_000, 20_000_000_000, 23_000_000_000)]
    native = dict(exit_code=0, timed_out=False, started_ns=1_000_000_000, finished_ns=22_000_000_000)
    measured = dict(native=native, points=points)
    context = dict(context=dict(observation_window=dict(scope_cpu_usec={"/user.slice": 10})))
    assert module.validate_trial(tmp_path, "idle", 123, 0, MANAGER, measured, context)["positive_control_response"] is None
    (tmp_path / "workload/unexpected").touch()
    with pytest.raises(ValueError, match="unexpected"):
        module.validate_trial(tmp_path, "idle", 123, 0, MANAGER, measured, context)


def test_fixed_order_and_owned_names():
    flat = [mode for block in ORDER for mode in block]
    assert len(flat) == 12 and all(flat.count(mode) == 3 for mode in set(flat))
    with pytest.raises(ValueError):
        unit_name(123, 12)
    with pytest.raises(ValueError):
        check_service_scope("0::/user.slice/other/owned.service", MANAGER, "owned.service")


def test_coordinator_starts_only_fixed_finite_user_service(tmp_path, monkeypatch):
    outcome = {}
    scope = MANAGER + "/app.slice/" + unit_name(123, 3)
    monkeypatch.setattr(module, "wait", lambda path, stopped: dict(cpus=list(range(20)))
        if path.name == "workload_ready.json" else dict(membership="0::" + scope + "\n"))
    commands = []

    def launch(command, **kwargs):
        commands.append(command)
        return SimpleNamespace(poll=lambda: 0, returncode=0)

    monkeypatch.setattr(module.subprocess, "Popen", launch)
    module.coordinate(tmp_path, "user-contended", 123, 3, MANAGER, threading.Event(), outcome)
    assert outcome["status"] == "completed"
    assert outcome["removal"][-1]["exists"] is False
    assert "--property=RuntimeMaxSec=45s" in commands[0]
    assert "--property=TasksMax=8" in commands[0]
    assert "--unit=orthohmm-root-123-3.service" in commands[0]
    assert json.loads((tmp_path / "workload_go.json").read_text()) == {"go": True}


@pytest.mark.parametrize("stop_code", [0, 1])
def test_failed_readiness_stops_only_its_owned_service(tmp_path, monkeypatch, stop_code):
    outcome, stopped_units = {}, []

    def wait(path, stopped):
        if path.name == "workload_ready.json":
            return dict(cpus=list(range(20)))
        raise RuntimeError("missing readiness")

    class Process:
        returncode = None
        def poll(self):
            return self.returncode
        def wait(self, timeout):
            self.returncode = -15
            return self.returncode

    monkeypatch.setattr(module, "wait", wait)
    monkeypatch.setattr(module.subprocess, "Popen", lambda *a, **k: Process())
    monkeypatch.setattr(module.subprocess, "run", lambda command, **k:
        stopped_units.append(command) or SimpleNamespace(returncode=stop_code, stdout="", stderr=""))
    module.coordinate(tmp_path, "user-contended", 123, 3, MANAGER, threading.Event(), outcome)
    assert outcome["status"] == "failed"
    assert stopped_units == [["systemctl", "--user", "stop", "orthohmm-root-123-3.service"]]
    assert not (tmp_path / "workload_go.json").exists()
    assert ("cleanup_error" in outcome) == bool(stop_code)

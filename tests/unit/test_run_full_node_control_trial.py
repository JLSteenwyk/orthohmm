import threading

import pytest

from benchmark_tools import run_full_node_control_trial as module


def test_cancelled_wait_terminates(tmp_path):
    stopped = threading.Event()
    stopped.set()
    with pytest.raises(RuntimeError, match="cancelled"):
        module.wait(tmp_path / "missing.json", stopped)


def test_readiness_timeout(tmp_path):
    with pytest.raises(TimeoutError, match="Missing control"):
        module.wait(tmp_path / "missing.json", threading.Event(), seconds=0)


@pytest.mark.parametrize("mode", ["steady", "churn"])
def test_controller_releases_native_workers(tmp_path, mode):
    module.save(tmp_path / "workload_ready.json", {"cpus": list(range(20))})
    outcome = {}
    module.coordinate(tmp_path, mode, threading.Event(), outcome)
    assert outcome == {"status": "completed"}
    assert module.read(tmp_path / "workload_go.json") == {"go": True}


def test_controller_cleanup_reaps_only_owned_competitor(tmp_path, monkeypatch):
    module.save(tmp_path / "workload_ready.json", {"cpus": list(range(20))})
    calls = []

    class Process:
        returncode = None

        def poll(self):
            return self.returncode

        def kill(self):
            calls.append("kill")
            self.returncode = -9

        def wait(self):
            calls.append("wait")
            return self.returncode

    monkeypatch.setattr(module.subprocess, "Popen", lambda *a, **kw: Process())
    original_wait = module.wait

    def fail(path, stopped, seconds=30):
        if path.name == "competitor_ready.json":
            raise TimeoutError("injected")
        return original_wait(path, stopped, seconds)

    monkeypatch.setattr(module, "wait", fail)
    outcome = {}
    module.coordinate(tmp_path, "contended", threading.Event(), outcome)
    assert outcome["status"] == "failed"
    assert outcome["competitor_exit_code"] == -9
    assert calls == ["kill", "wait"]
    assert not (tmp_path / "workload_go.json").exists()


def test_collector_failure_cancels_and_joins_controller(tmp_path, monkeypatch):
    ended = threading.Event()

    def coordinate(directory, mode, stopped, outcome):
        stopped.wait()
        outcome["status"] = "cancelled"
        ended.set()

    def measure(*args):
        raise ValueError("injected collector failure")

    monkeypatch.setattr(module, "coordinate", coordinate)
    monkeypatch.setattr(module, "measure", measure)
    with pytest.raises(ValueError, match="collector failure"):
        module.trial(tmp_path / "trial", "steady", 1)
    assert ended.is_set()
    assert module.read(tmp_path / "trial/controller.json") == {"status": "cancelled"}


def test_invalid_mode_does_not_create_directory(tmp_path):
    path = tmp_path / "trial"
    with pytest.raises(ValueError, match="Unknown"):
        module.trial(path, "invalid", 1)
    assert not path.exists()


@pytest.mark.parametrize("mode,detected", [("steady", False), ("churn", True),
                                         ("contended", True), ("contended", False)])
def test_complete_trial_wiring(tmp_path, monkeypatch, mode, detected):
    # Synthetic collector/worker evidence tests orchestration, not CPU validity.
    screens = [dict(screen_passed=not detected,
                    reasons=["excess_unassigned_cpu"] if detected else [])]
    report = dict(points=[dict(native_membership="native", host=[dict(raw=dict(cgroup_membership="batch"))])],
                  native={"fixture": "native"}, screening=dict(narrow_intervals=screens))
    validation = dict(status="fixture_validation", common_started_ns=1, common_finished_ns=2)
    calls = []

    def coordinate(directory, condition, stopped, outcome):
        module.save(directory / "workload_go.json", {"go": True})
        for name in ("workload_ready", "workload_done"):
            module.save(directory / f"{name}.json", {"fixture": name})
        if condition == "contended":
            for name in ("competitor_ready", "competitor_done"):
                module.save(directory / f"{name}.json", {"fixture": name})
            outcome["competitor_exit_code"] = 0
        outcome["status"] = "completed"

    def measure(command, directory, job, cpus, memory, timeout, interval):
        assert command[-2:] == ["--worker", "churn" if mode == "churn" else "steady"]
        assert (job, cpus, memory, timeout, interval) == (17, 20, 96*1024**3, 60, 1.)
        return report

    def validate(*args):
        calls.append(args)
        return validation

    monkeypatch.setattr(module, "coordinate", coordinate)
    monkeypatch.setattr(module, "measure", measure)
    monkeypatch.setattr(module, "evaluate", lambda *args: report["screening"])
    monkeypatch.setattr(module, "validate_witnesses", validate)
    monkeypatch.setattr(module, "common_intervals", lambda points, witness: [0])
    result = module.trial(tmp_path / "trial", mode, 17)
    assert result["positive_control_detected"] == (detected if mode == "contended" else None)
    assert result["common_narrow_flagged"] == ([0] if detected else [])
    assert result["status"] == "workload_validated"
    assert result["scientific_timings_admitted"] is False
    assert len(calls) == 1
    assert calls[0][0] == mode
    assert calls[0][4:6] == ("native", "batch")
    assert module.read(tmp_path / "trial/trial.json") == result

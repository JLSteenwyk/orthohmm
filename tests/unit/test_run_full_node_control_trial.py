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

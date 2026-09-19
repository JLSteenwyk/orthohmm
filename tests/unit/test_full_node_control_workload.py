import json
import os

import pytest

from benchmark_tools import full_node_control_workload as module


@pytest.mark.parametrize("mode", ["steady", "churn"])
def test_bounded_two_worker_local_smoke(tmp_path, mode):
    cpus = sorted(os.sched_getaffinity(0))[:2]
    module.save(tmp_path / "workload_go.json", {"go": True})
    result = module.workload(tmp_path, mode, cpus, duration=.1)
    assert set(result["statuses"].values()) == {0}
    assert len(result["workers"]) == len(cpus)
    for row, cpu in zip(result["workers"], cpus):
        assert row["affinity"] == row["final_affinity"] == [cpu]
        assert row["membership"] == row["final_membership"] == result["membership"]
        assert row["parent_pid"] == os.getpid()
        assert row["finished_ns"] - row["started_ns"] >= 100000000
        assert not row["creation_cap_reached"]
        assert row["creations"] > 0 if mode == "churn" else row["creations"] == 0
        with pytest.raises(ChildProcessError):
            os.waitpid(row["pid"], os.WNOHANG)
    assert json.loads((tmp_path / "workload_done.json").read_text()) == result
    assert not result["scientific_timings_admitted"]
    with pytest.raises(FileExistsError):
        module.workload(tmp_path, mode, cpus, duration=.1)


def test_creation_cap_preserves_invalid_workload_evidence(tmp_path):
    cpu = min(os.sched_getaffinity(0))
    module.save(tmp_path / "workload_go.json", {"go": True})
    with pytest.raises(RuntimeError, match="cap reached"):
        module.workload(tmp_path, "churn", [cpu], duration=.1, cap=1)
    result = json.loads((tmp_path / "workload_done.json").read_text())
    assert result["workers"][0]["creations"] == 1
    assert result["workers"][0]["creation_cap_reached"]


def test_failed_worker_is_reaped(tmp_path):
    cpu = min(os.sched_getaffinity(0))
    module.save(tmp_path / "workload_go.json", {"go": False})
    with pytest.raises(RuntimeError, match="worker failed"):
        module.workload(tmp_path, "steady", [cpu], duration=.1)
    ready = json.loads((tmp_path / f"ready_{cpu}.json").read_text())
    assert (tmp_path / f"failed_{cpu}.json").exists()
    with pytest.raises(ChildProcessError):
        os.waitpid(ready["pid"], os.WNOHANG)
    assert not (tmp_path / "workload_done.json").exists()


@pytest.mark.parametrize("field,value", [("mode", "invalid"), ("cpus", []), ("cpus", [True]),
    ("duration", 0), ("duration", 21), ("duration", float("nan")), ("cap", 0), ("cap", 200001)])
def test_invalid_resource_bounds_rejected(field, value):
    args = dict(mode="steady", cpus=[min(os.sched_getaffinity(0))], duration=.1, cap=10)
    args[field] = value
    with pytest.raises(ValueError):
        module.validate(**args)


def test_lost_parent_stops_work():
    with pytest.raises(RuntimeError, match="parent disappeared"):
        module.work("steady", .1, 10, -1)


def test_competitor_local_smoke_restores_test_affinity(tmp_path):
    allowed = os.sched_getaffinity(0)
    cpu = min(allowed)
    module.save(tmp_path / "workload_go.json", {"go": True})
    try:
        result = module.competitor(tmp_path, cpu, duration=.02)
        assert result["affinity"] == result["final_affinity"] == [cpu]
        assert result["creations"] == 0
        assert result["self_cpu_s"] > 0
    finally:
        os.sched_setaffinity(0, allowed)


@pytest.mark.parametrize("value", [{"go": 1}, {"go": False}, {"go": True, "extra": 0}, [], None])
def test_start_signal_requires_exact_boolean(value):
    with pytest.raises(ValueError, match="shared start"):
        module.require_start(value)


@pytest.mark.parametrize("failure", ["fork", "readiness"])
def test_partial_startup_reaps_owned_children(tmp_path, monkeypatch, failure):
    cpus = sorted(os.sched_getaffinity(0))[:2]
    if len(cpus) < 2:
        pytest.skip("Needs two allowed CPUs for partial startup")
    original_fork = os.fork
    original_wait = module.wait_file
    owned = []

    def tracked_fork():
        if failure == "fork" and owned:
            raise OSError("injected fork failure")
        pid = original_fork()
        if pid:
            owned.append(pid)
        return pid

    def failed_readiness(path, seconds):
        if failure == "readiness" and path.name.startswith("ready_"):
            raise TimeoutError("injected readiness failure")
        return original_wait(path, seconds=seconds)

    monkeypatch.setattr(module.os, "fork", tracked_fork)
    monkeypatch.setattr(module, "wait_file", failed_readiness)
    with pytest.raises((OSError, TimeoutError), match="injected"):
        module.workload(tmp_path, "steady", cpus, duration=.1)
    assert len(owned) == (1 if failure == "fork" else 2)
    for pid in owned:
        with pytest.raises(ChildProcessError):
            os.waitpid(pid, os.WNOHANG)
    assert not (tmp_path / "workload_done.json").exists()

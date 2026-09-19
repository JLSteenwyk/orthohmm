import os
import subprocess

import pytest

from benchmark_tools import cpu_creation_workload as module


@pytest.mark.parametrize("mode", ["steady", "creation"])
def test_bounded_real_worker_and_child_accounting(mode):
    before = os.sched_getaffinity(0)
    cpu = min(before)
    result = module.run(mode, [cpu], .08, .001, 100)
    assert os.sched_getaffinity(0) == before
    assert result["scientific_timings_admitted"] is False
    row = result["workers"][0]
    assert row["affinity"] == [cpu] and row["finished_ns"] > row["started_ns"]
    assert result["worker_cpu_s"] > 0
    if mode == "creation":
        assert 0 < row["reaped_children"] <= 100
        assert row["children_cpu_s"] > 0
    else:
        assert row["reaped_children"] == 0 and row["children_cpu_s"] == 0


def test_child_count_cap_is_retained():
    result = module.run("creation", [min(os.sched_getaffinity(0))], .3, .001, 1)
    assert result["reaped_children"] == 1
    assert result["child_count_cap_reached"] is True


@pytest.mark.parametrize("mode,seconds,child,cap", [
    ("unknown", 1, .005, 10), ("steady", 0, .005, 10), ("steady", 31, .005, 10),
    ("steady", float("nan"), .005, 10), ("creation", 1, 0, 10),
    ("creation", 1, float("inf"), 10), ("creation", 1, .005, 10001),
    ("creation", 1, .005, True)])
def test_limits_reject_before_launch(mode, seconds, child, cap):
    with pytest.raises(ValueError):
        module.run(mode, [min(os.sched_getaffinity(0))], seconds, child, cap)


@pytest.mark.parametrize("cpus", [[], [-1], [True], list(range(21))])
def test_invalid_affinity(cpus):
    with pytest.raises(ValueError):
        module.run("steady", cpus, .01)


def test_duplicate_cpu():
    cpu = min(os.sched_getaffinity(0))
    with pytest.raises(ValueError):
        module.run("steady", [cpu, cpu], .01)


def test_timeout_reaps_launched_worker(monkeypatch):
    real_popen = subprocess.Popen
    processes = []

    def launch(*args, **kwargs):
        process = real_popen(*args, **kwargs)
        processes.append(process)
        def timeout(**unused):
            raise subprocess.TimeoutExpired(args[0], .001)
        monkeypatch.setattr(process, "communicate", timeout)
        return process

    monkeypatch.setattr(module.subprocess, "Popen", launch)
    with pytest.raises(subprocess.TimeoutExpired):
        module.run("creation", [min(os.sched_getaffinity(0))], .1)
    assert len(processes) == 1 and processes[0].poll() is not None
    assert processes[0].stdout.closed and processes[0].stderr.closed

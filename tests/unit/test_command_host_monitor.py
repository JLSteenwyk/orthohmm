import io
import json

import pytest

from benchmark_tools.command_host_monitor import HostMonitor


def sample(t, cpu=0., errors=None):
    return {"started_monotonic_s": t, "finished_monotonic_s": t + .01, "errors": errors or [],
            "processes": [{"pid": 999999, "created": 1., "cgroup": "/other", "name": "fixture",
                           "user_s": cpu, "system_s": 0., "observed_monotonic_s": t}]}


@pytest.mark.parametrize("cpu,errors,state", [(1., [], "competing_cpu_observed"),
    (0., [], "no_large_persistent_competitor_observed"), (0., [{}], "inconclusive")])
def test_bracketed_stream_retains_intervals(cpu, errors, state):
    sequence = iter([sample(1.), sample(3., cpu, errors)])
    handle = io.StringIO()
    monitor = HostMonitor(handle, "/job/target", lambda: next(sequence))
    monitor.observe()
    monitor.observe()
    result = monitor.summary(1.5, 2.5)
    assert result["status"] == state
    assert result["command_bracketed_by_samples"]
    assert result["controlled_workload_verified"] is False
    assert len(handle.getvalue().splitlines()) == 2


def test_observation_failure_cannot_become_quiet():
    sequence = iter([sample(1.), None, sample(3.), sample(4.)])
    def next_sample():
        value = next(sequence)
        if value is None:
            raise OSError("fixture")
        return value
    handle = io.StringIO()
    monitor = HostMonitor(handle, "/target", next_sample)
    for _ in range(4):
        monitor.observe()
    result = monitor.summary(1.5, 3.5)
    assert result["status"] == "inconclusive" and result["observation_errors"] == 1
    assert result["command_bracketed_by_samples"]
    assert json.loads(handle.getvalue().splitlines()[1])["observation_error"] == "OSError"


@pytest.mark.parametrize("launch,end", [(0., 2.5), (1.5, 4.)])
def test_partial_coverage_inconclusive(launch, end):
    sequence = iter([sample(1.), sample(3.)])
    monitor = HostMonitor(io.StringIO(), "/target", lambda: next(sequence))
    monitor.observe()
    monitor.observe()
    assert monitor.summary(launch, end)["status"] == "inconclusive"


def test_no_samples_safe_and_inconclusive():
    monitor = HostMonitor(io.StringIO(), "/target")
    assert monitor.summary(1., 2.)["status"] == "inconclusive"


def test_short_lived_process_churn_is_not_quiet():
    first, last = sample(1.), sample(3.)
    last["processes"] = []
    sequence = iter([first, last])
    monitor = HostMonitor(io.StringIO(), "/target", lambda: next(sequence))
    monitor.observe()
    monitor.observe()
    assert monitor.summary(1.5, 2.5)["status"] == "inconclusive"

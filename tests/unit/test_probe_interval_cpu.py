from copy import deepcopy

import pytest

from benchmark_tools.probe_host_counters import parse_cpu
from benchmark_tools.probe_interval_cpu import evaluate, interval, validate_point


def sample(start, busy):
    raw = f"cpu {busy} 0 0 1000 0 0 0 0 0 0\n"
    return dict(started_monotonic_ns=start, finished_monotonic_ns=start + 100,
        raw=dict(proc_stat=raw, cgroup_membership="0::/job_123/step_batch\n",
                 boot_id="boot", online_cpus="0-19"),
        cpu_ticks=parse_cpu(raw), errors=[], optional={})


def point(index, extra=0):
    start = index * 1_000_000_000
    return dict(host=[sample(start, index * 100 + extra), sample(start + 400, index * 100 + extra)],
        native_membership="0::/job_123/step_0/task_0\n", native_cpu_scope="/job_123/step_0",
        native_read_ns=[start + 200, start + 300], native_cpu_stat=f"usage_usec {index * 1000000}\n",
        ticks_per_second=100)


def test_whole_window_dilutes_burst_but_interval_flags_it():
    points = [point(i, extra=75 if i >= 3 else 0) for i in range(11)]
    result = evaluate(points, 123)
    assert result["whole_window"]["screen_passed"] is True
    assert result["whole_window"]["signed_unassigned_average_cores"] == .075
    assert result["flagged_intervals"] == [2]
    assert result["intervals"][2]["signed_unassigned_average_cores"] == .75
    assert result["controlled_workload_verified"] is False
    assert result["scientific_timings_admitted"] is False


def test_quiet_and_overhang():
    result = evaluate([point(i) for i in range(11)], 123)
    assert result["flagged_intervals"] == []
    assert result["intervals"][0]["outer_read_overhang_s"] == 4e-7


@pytest.mark.parametrize("fault", ["error", "scope", "read", "tick", "boot", "counter", "identity", "gap", "overlap"])
def test_invalid_interval_evidence(fault):
    left, right = point(0), point(1)
    if fault == "error":
        right["host"][1]["errors"] = ["missing"]
    elif fault == "scope":
        right["native_cpu_scope"] += "/task_0"
    elif fault == "read":
        right["native_read_ns"][0] -= 500
    elif fault == "tick":
        right["ticks_per_second"] = 0
    elif fault == "boot":
        right["host"][1]["raw"]["boot_id"] = "other"
    elif fault == "counter":
        left["native_cpu_stat"] = "usage_usec 2000000\n"
    elif fault == "identity":
        right["native_membership"] = "0::/job_123/step_1/task_0\n"
        right["native_cpu_scope"] = "/job_123/step_1"
    elif fault == "gap":
        right = point(2)
    else:
        right = deepcopy(left)
    with pytest.raises(ValueError):
        interval(left, right, 123)


def test_missing_observation_is_not_interpolated():
    with pytest.raises(ValueError, match="Incomplete"):
        evaluate([point(i) for i in range(10)], 123)


def test_wrong_job_rejected():
    with pytest.raises(ValueError, match="job scope"):
        validate_point(point(0), 999)


def test_negative_accounting_preserved():
    right = point(1)
    right["native_cpu_stat"] = "usage_usec 2000000\n"
    result = interval(point(0), right, 123)
    assert result["signed_unassigned_cpu_s"] == -1
    assert result["reasons"] == ["negative_accounting_discrepancy"]

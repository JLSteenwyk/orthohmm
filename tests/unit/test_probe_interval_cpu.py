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


def test_actual_dgx_control_replay_and_burst_boundaries():
    import hashlib
    import json
    from pathlib import Path

    root = Path(__file__).resolve().parents[2]
    evidence = root / "benchmark_tools/results/dgx_interval_controls_21806.json"
    report = json.loads(evidence.read_text())
    assert report["job_id"] == 21806
    assert report["all_control_expectations_met"] is True
    assert report["scientific_timings_admitted"] is False
    assert report["controlled_workload_verified"] is False
    for name, expected in report["sources"].items():
        assert hashlib.sha256((root / "benchmark_tools" / name).read_bytes()).hexdigest() == expected
    assert [t["mode"] for t in report["trials"]] == ["quiet", "completed_burst"]
    for trial in report["trials"]:
        points = trial["points"]
        assert evaluate(points, 21806) == trial["result"]
        assert trial["native"]["started_ns"] > points[0]["host"][1]["finished_monotonic_ns"]
        assert trial["native"]["finished_ns"] > points[-1]["host"][1]["finished_monotonic_ns"]
        assert trial["native"]["cpu_s"] > 9
        assert trial["result"]["whole_window"]["screen_passed"] is True
    quiet, burst = report["trials"]
    assert quiet["result"]["flagged_intervals"] == []
    assert burst["result"]["flagged_intervals"] == [2]
    assert burst["result"]["intervals"][2]["reasons"] == ["excess_unassigned_cpu"]
    load = burst["sibling"]
    assert .75 <= load["cpu_seconds"] < .85
    assert burst["points"][2]["host"][1]["finished_monotonic_ns"] < load["started_ns"]
    assert load["started_ns"] < load["finished_ns"] < burst["points"][3]["host"][0]["started_monotonic_ns"]
    assert load["cgroup"] == burst["points"][0]["host"][0]["raw"]["cgroup_membership"]
    scheduler = (root / "benchmark_tools/results/dgx_interval_controls_21806_scheduler.txt").read_text().split()
    for field in ("JobId=21806", "JobState=COMPLETED", "ExitCode=0:0", "Restarts=0",
                  "NodeList=spark-7ff0", "NumCPUs=20", "CPUs/Task=2", "OverSubscribe=NO"):
        assert field in scheduler

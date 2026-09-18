import copy

import pytest

from benchmark_tools.probe_dgx_cpu_hierarchy import compare, usage, validate_point
from benchmark_tools.probe_host_counters import parse_cpu


def point(i):
    start = i * 1_000_000_000
    def host(offset):
        raw = f"cpu {100*i} 0 0 1000 0 0 0 0 0 0\n"
        return dict(started_monotonic_ns=start+offset, finished_monotonic_ns=start+offset+10,
                    raw=dict(proc_stat=raw, cgroup_membership="0::/job_123/step_batch/task_0\n",
                             boot_id="boot", online_cpus="0-19"), cpu_ticks=parse_cpu(raw), errors=[])
    def cpu(scope, offset, value):
        return dict(scope=scope, started_ns=start+offset, finished_ns=start+offset+5,
                    raw=f"usage_usec {value}\n")
    return dict(host=[host(0), host(80)], ticks=100, native_membership="0::/job_123/step_0/task_0\n",
                parent=[cpu("/job_123", 20, i*1000000), cpu("/job_123", 60, i*1000000+100)],
                children=[cpu("/job_123/step_0", 30, i*100000), cpu("/job_123/step_batch", 40, i*750000)])


def test_disjoint_steps_and_signed_residuals():
    result = compare(point(0), point(1), 123)
    assert result["step_cpu_s"] == dict(step_0=.1, step_batch=.75)
    assert result["job_outer_cpu_s"] == 1.0001
    assert result["job_inner_cpu_s"] == .9999
    assert result["job_outer_minus_step_sum_cpu_s"] == pytest.approx(.1501)
    assert result["host_minus_job_outer_cpu_s"] == pytest.approx(-.0001)
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["nested", "duplicate", "wrong_job", "missing_batch", "order", "enclosure", "error", "parent", "decrease"])
def test_bad_point_rejected(fault):
    p = point(0)
    if fault == "nested":
        p["children"][0]["scope"] += "/task_0"
    elif fault == "duplicate":
        p["children"].append(copy.deepcopy(p["children"][0]))
    elif fault == "wrong_job":
        p["native_membership"] = "0::/job_999/step_0/task_0\n"
    elif fault == "missing_batch":
        p["children"].pop()
    elif fault == "order":
        p["children"].reverse()
    elif fault == "enclosure":
        p["parent"][0]["started_ns"] = 0
    elif fault == "error":
        p["host"][0]["errors"] = ["failed"]
    elif fault == "parent":
        p["parent"][0]["scope"] = "/"
    else:
        p["parent"][0]["raw"] = "usage_usec 101\n"
    with pytest.raises(ValueError):
        validate_point(p, 123)


@pytest.mark.parametrize("fault", ["topology", "child_decrease", "overlap", "membership"])
def test_bad_interval_rejected(fault):
    a, b = point(0), point(1)
    if fault == "topology":
        b["children"].append(dict(scope="/job_123/step_extern", started_ns=1000000050,
                                  finished_ns=1000000055, raw="usage_usec 0\n"))
    elif fault == "child_decrease":
        a["children"][0]["raw"] = "usage_usec 100001\n"
    elif fault == "membership":
        b["native_membership"] = "0::/job_123/step_0/task_1\n"
    else:
        b = copy.deepcopy(a)
    with pytest.raises(ValueError):
        compare(a, b, 123)


@pytest.mark.parametrize("text", ["usage_usec -1\n", "usage_usec 1\nusage_usec 2\n", "usage_usec nan\n"])
def test_malformed_raw_counter_rejected(text):
    with pytest.raises(ValueError):
        usage(text)


def test_actual_dgx_control_replay_and_source_identity():
    import hashlib
    import json
    from pathlib import Path

    root = Path(__file__).resolve().parents[2]
    report = json.loads((root / "benchmark_tools/results/dgx_cpu_hierarchy_controls_21816.json").read_text())
    assert report["job_id"] == 21816
    assert report["scientific_timings_admitted"] is False
    assert report["controlled_workload_verified"] is False
    assert [t["mode"] for t in report["trials"]] == ["quiet", "completed_burst", "sustained_batch"]
    for name, sha in report["sources"].items():
        assert hashlib.sha256((root / "benchmark_tools" / name).read_bytes()).hexdigest() == sha
    for trial, repeats in zip(report["trials"], (0, 1, 4)):
        assert compare(*trial["points"], 21816) == trial["result"]
        assert len(trial["loads"]) == repeats
        if repeats:
            assert trial["load_localization_expectation_met"] is True
            assert trial["result"]["step_cpu_s"]["step_batch"] >= .5 * repeats
        else:
            assert trial["load_localization_expectation_met"] is None
        for load in trial["loads"]:
            assert .75 <= load["cpu_seconds"] < .85
            assert load["cgroup"] == trial["points"][0]["host"][0]["raw"]["cgroup_membership"]
    scheduler = (root / "benchmark_tools/results/dgx_cpu_hierarchy_controls_21816_scheduler.txt").read_text().split()
    for token in ("JobState=COMPLETED", "ExitCode=0:0", "Restarts=0", "NodeList=spark-7ff0",
                  "NumCPUs=20", "CPUs/Task=2", "OverSubscribe=NO", "MinMemoryNode=2G"):
        assert token in scheduler

from copy import deepcopy

import pytest

from benchmark_tools.replay_dgx_resource_samples import validate_rows
from benchmark_tools.monitor_slurm_resources import summarize


def fixture():
    scope = "/job_123/step_batch/task"
    rows = []
    for index, begin in enumerate((1., 3., 5.)):
        cpu = {"usage_usec": 30 + index * 30, "user_usec": 20 + index * 20, "system_usec": 10 + index * 10}
        metrics = {"cpu": cpu, "memory_events": {"oom": 0}, "memory_stat": {"anon": 10},
            "effective_cpus": [0, 1], "memory_current_bytes": 10,
            "memory_peak_since_creation_or_reset_bytes": 20, "memory_max_bytes": 100}
        metrics["raw"] = {"cpu.stat": "\n".join(f"{k} {v}" for k, v in cpu.items()),
            "memory.events": "oom 0", "memory.stat": "anon 10", "cpuset.cpus.effective": "0-1",
            "memory.current": "10", "memory.peak": "20", "memory.max": "100"}
        rows.append({"index": index, "anchor_created": 100., "elapsed_s": begin + .1 - .5,
            "started_monotonic_s": begin, "finished_monotonic_s": begin + .1,
            "snapshot": {"job_id": 123, "pid": 9, "scope": scope, "proc_cgroup": "0::" + scope,
                "status": "resource_snapshot", "metrics": metrics, "sampling_errors": [],
                "sampled_sum_process_rss_bytes": 12, "ancestor_limits": [{"memory_max": "100"}],
                "processes": [{"pid": 9, "created": 100., "rss_bytes": 12, "cgroup": scope}]}})
    measured = {"status": "command_exited_zero", "timed_out": False,
        "clock_domain": {"hostname": "spark-7ff0", "clock": "time.monotonic", "unit": "seconds", "boot_id": "fixture"},
        "wrapper_started_monotonic_s": .5, "command_launch_started_monotonic_s": 2.,
        "command_wait_finished_monotonic_s": 4., "wrapper_finished_monotonic_s": 6.,
        "command_wall_s": 2., "wrapper_wall_s": 5.5, "wrapper_pid": 9, "job_id": 123,
        "baseline_scope": scope, "requested_cpus": 2, "requested_memory_bytes": 100,
        "summary": summarize(rows)}
    return measured, rows


def test_reproduces_exact_resource_summary():
    measured, rows = fixture()
    assert validate_rows(measured, rows) == measured["summary"]


@pytest.mark.parametrize("problem", ["clock", "wall", "index", "order", "elapsed", "identity", "raw", "rss", "duplicate",
    "cpu", "bracket", "leftover", "summary", "counter_decrease"])
def test_corruption_fails(problem):
    measured, rows = fixture()
    if problem == "clock":
        measured["clock_domain"]["clock"] = "time.time"
    elif problem == "wall":
        measured["command_wall_s"] = 1.
    elif problem == "index":
        rows[1]["index"] = 9
    elif problem == "order":
        rows[1]["started_monotonic_s"] = .1
    elif problem == "elapsed":
        rows[1]["elapsed_s"] += 1
    elif problem == "identity":
        rows[1]["snapshot"]["job_id"] = 124
    elif problem == "raw":
        rows[1]["snapshot"]["metrics"]["raw"]["memory.current"] = "99"
    elif problem == "rss":
        rows[1]["snapshot"]["sampled_sum_process_rss_bytes"] = 100
    elif problem == "duplicate":
        rows[1]["snapshot"]["processes"] *= 2
    elif problem == "cpu":
        measured["requested_cpus"] = 3
    elif problem == "bracket":
        measured["command_launch_started_monotonic_s"] = .9
        measured["command_wall_s"] = 3.1
    elif problem == "leftover":
        rows[-1]["snapshot"]["processes"][0]["pid"] = 10
    elif problem == "summary":
        measured["summary"]["maximum_sampled_sum_rss_bytes"] = 99
    elif problem == "counter_decrease":
        rows[-1]["snapshot"]["metrics"] = deepcopy(rows[0]["snapshot"]["metrics"])
    with pytest.raises(ValueError):
        validate_rows(measured, rows)

import json
import os
import subprocess
import sys
import time

import psutil

import pytest

from benchmark_tools.measure_slurm_command import check_allocation, measure, stop_owned_group
from benchmark_tools import measure_slurm_command as measurement


def snapshot(pid, job):
    return {"pid": pid, "job_id": job, "scope": "/job_123/step_batch/user/task_0",
            "processes": [{"pid": pid, "rss_bytes": 100}], "sampling_errors": [], "sampled_sum_process_rss_bytes": 100,
            "metrics": {"effective_cpus": [0], "cpu": {"usage_usec": 1, "user_usec": 1, "system_usec": 0},
                        "memory_peak_since_creation_or_reset_bytes": 100},
            "ancestor_limits": [{"memory_max": "1024"}]}


@pytest.mark.parametrize("problem", [None, "process", "errors", "cpus", "memory", "unlimited"])
def test_allocation_gate(problem):
    sample = snapshot(os.getpid(), 123)
    if problem == "process":
        sample["processes"].append({"pid": os.getpid() + 1})
    elif problem == "errors":
        sample["sampling_errors"] = [{}]
    elif problem == "cpus":
        sample["metrics"]["effective_cpus"] = [0, 1]
    elif problem == "memory":
        sample["ancestor_limits"] = [{"memory_max": "2048"}]
    elif problem == "unlimited":
        sample["ancestor_limits"] = [{"memory_max": "max"}]
    if problem:
        with pytest.raises(ValueError):
            check_allocation(sample, os.getpid(), 1, 1024)
    else:
        check_allocation(sample, os.getpid(), 1, 1024)


@pytest.mark.parametrize("code,expected", [(0, "command_exited_zero"), (7, "command_failed")])
def test_real_command_exit_and_retained_logs(tmp_path, code, expected):
    result = measure([sys.executable, "-c", f"print('fixture'); raise SystemExit({code})"], tmp_path / "run", 123, 1, 1024, 5., .02, snapshot)
    assert result["status"] == expected and result["exit_code"] == code
    assert result["summary"]["observations"] >= 2
    assert "fixture" in (tmp_path / "run/command.log").read_text()
    rows = [json.loads(line) for line in (tmp_path / "run/samples.jsonl").read_text().splitlines()]
    start, end = result["command_launch_started_monotonic_s"], result["command_wait_finished_monotonic_s"]
    assert result["wrapper_started_monotonic_s"] <= rows[0]["started_monotonic_s"]
    assert rows[0]["finished_monotonic_s"] <= start <= end <= rows[-1]["started_monotonic_s"]
    assert rows[-1]["finished_monotonic_s"] <= result["wrapper_finished_monotonic_s"]
    assert all(row["started_monotonic_s"] <= row["finished_monotonic_s"] for row in rows)
    assert result["command_wall_s"] == end - start
    assert result["clock_domain"]["boot_id"]


def test_timeout_stops_owned_command(tmp_path):
    result = measure([sys.executable, "-c", "import time; time.sleep(30)"], tmp_path / "timeout", 123, 1, 1024, .1, .02, snapshot)
    assert result["status"] == "command_timed_out" and result["timed_out"]
    assert result["exit_code"] < 0


def test_sampling_failure_retains_native_success(tmp_path):
    calls = 0
    def failing(pid, job):
        nonlocal calls
        calls += 1
        if calls > 1:
            raise OSError("fixture sampling failure")
        return snapshot(pid, job)
    result = measure([sys.executable, "-c", "import time; time.sleep(.15); print('finished')"], tmp_path / "run", 123, 1, 1024, 5., .02, failing)
    assert result["status"] == "measurement_failed" and result["exit_code"] == 0
    assert "finished" in (tmp_path / "run/command.log").read_text()


def test_spawn_failure_is_recorded(tmp_path):
    with pytest.raises(FileNotFoundError):
        measure([str(tmp_path / "absent")], tmp_path / "run", 123, 1, 1024, 5., .02, snapshot)
    assert json.loads((tmp_path / "run/results.json").read_text())["status"] == "wrapper_failed"


def test_refuse_existing_output(tmp_path):
    with pytest.raises(FileExistsError):
        measure(["unused"], tmp_path, 123, 1, 1024, 5., snapshot_fn=snapshot)


@pytest.mark.parametrize("failure", [False, True])
def test_host_monitor_does_not_change_native_completion(tmp_path, failure):
    def host_snapshot():
        if failure:
            raise OSError("host fixture failure")
        now = time.monotonic()
        return {"started_monotonic_s": now, "finished_monotonic_s": now,
                "processes": [], "errors": []}
    result = measure([sys.executable, "-c", "import time; time.sleep(.1); print('finished')"],
                     tmp_path / "run", 123, 1, 1024, 5., .02, snapshot,
                     monitor_host=True, host_snapshot_fn=host_snapshot)
    assert result["status"] == "command_exited_zero" and result["exit_code"] == 0
    assert result["host_workload"]["status"] == ("inconclusive" if failure else "no_large_persistent_competitor_observed")
    assert result["controlled_workload_verified"] is False
    assert result["host_workload"]["command_bracketed_by_samples"] is not failure
    assert "host_samples.jsonl" in result
    assert result["host_workload"]["observer_pid"] == result["wrapper_pid"] == os.getpid()
    for key in ("command_launch_started_monotonic_s", "command_wait_finished_monotonic_s"):
        assert result["host_workload"][key] == result[key]


def test_timeout_cleanup_reaches_child_that_ignores_term():
    child_code = "import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); print('ready',flush=True); time.sleep(30)"
    parent_code = ("import subprocess,sys,time; "
                   f"p=subprocess.Popen([sys.executable,'-c',{child_code!r}],stdout=subprocess.PIPE,text=True); "
                   "p.stdout.readline(); print(p.pid,flush=True); time.sleep(30)")
    process = subprocess.Popen([sys.executable, "-c", parent_code], start_new_session=True, stdout=subprocess.PIPE, text=True)
    try:
        child_pid = int(process.stdout.readline())
        stop_owned_group(process, grace=.05)
        for _ in range(100):
            try:
                child = psutil.Process(child_pid)
                if child.status() == psutil.STATUS_ZOMBIE:
                    break
            except psutil.NoSuchProcess:
                break
            time.sleep(.01)
        else:
            pytest.fail("Owned child remained active after timeout cleanup")
    finally:
        stop_owned_group(process, grace=.05)
        process.stdout.close()


@pytest.mark.parametrize("interval,expected", [(30., 2), (1., 5)])
def test_host_cadence_independent_of_resource_sampling(tmp_path, monkeypatch, interval, expected):
    clock = [0.]

    class Command:
        pid = 999999
        returncode = None
        waits = 0

        def wait(self, timeout):
            self.waits += 1
            clock[0] += 1.
            if self.waits < 4:
                raise subprocess.TimeoutExpired("fixture", timeout)
            self.returncode = 0
            return 0

    def host_snapshot():
        return {"started_monotonic_s": clock[0], "finished_monotonic_s": clock[0],
                "processes": [], "errors": []}

    monkeypatch.setattr(measurement.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(measurement.subprocess, "Popen", lambda *a, **kw: Command())
    result = measure(["fixture"], tmp_path / "run", 123, 1, 1024, 10., .02, snapshot,
                     monitor_host=True, host_snapshot_fn=host_snapshot, host_interval_s=interval)
    assert result["status"] == "command_exited_zero"
    assert result["summary"]["observations"] == 5
    assert result["host_interval_s"] == interval
    assert result["host_workload"]["successful_snapshots"] == expected
    assert result["host_workload"]["command_bracketed_by_samples"]


@pytest.mark.parametrize("interval", [0., -1., float("inf"), float("nan")])
def test_invalid_host_cadence_rejected_before_launch(tmp_path, interval):
    with pytest.raises(ValueError):
        measure(["unused"], tmp_path / "absent", 123, 1, 1024, 5., host_interval_s=interval)
    assert not (tmp_path / "absent").exists()

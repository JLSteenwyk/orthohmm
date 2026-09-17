import json
import os
import subprocess
import sys
import time

import psutil

import pytest

from benchmark_tools.measure_slurm_command import check_allocation, measure, stop_owned_group


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

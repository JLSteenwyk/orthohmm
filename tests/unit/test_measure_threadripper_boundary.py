import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import measure_threadripper_boundary as module


def setup(tmp_path, monkeypatch, *, cycles=1, code=0, timed_out=False, fault=None):
    directory = tmp_path / "measurement"
    events, points = [], []
    clock = [0.]
    for key, value in dict(SLURM_JOB_ID="42", SLURM_CPUS_PER_TASK="64", SLURM_MEM_PER_NODE="131072").items():
        monkeypatch.setenv(key, value)
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    monkeypatch.setattr(module.time, "monotonic", lambda: clock[0])
    wall = 85800. if timed_out else float(cycles)
    done = dict(exit_code=code, timed_out=timed_out, started_ns=100,
                finished_ns=100 + int(wall * 1e9))

    class Process:
        finished = False
        def poll(self):
            return 0 if self.finished or fault == "early_exit" and clock[0] > 0 else None
        def wait(self, timeout):
            assert (directory / "release.json").exists()
            self.finished = True
            events.append("cleanup")
            return 7 if fault == "wrapper" else 0

    def launch(argv, **kwargs):
        assert argv[:7] == ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1",
            "--cpus-per-task=64", "--cpu-bind=mask_cpu:0xffffffff"]
        assert argv[-3:] == [str(Path(module.periodic.__file__).resolve()), "--worker", str(directory)]
        return Process()

    monkeypatch.setattr(module.subprocess, "Popen", launch)
    monkeypatch.setattr(module, "wait_file", lambda *a: dict(pid=1,
        cgroup="0::/slurm/job_42/step_0/user/task_0\n", placement={}))
    host = SimpleNamespace(observe=lambda: events.append("host_initial"))
    monkeypatch.setattr(module, "HostMonitor", lambda *a, **k: host)
    class Observer:
        def __init__(self, actual_host, period):
            assert actual_host is host and period == 30.
        def start(self):
            assert not (directory / "go.json").exists()
            events.append("host_start")
        def finish(self, start, end):
            assert start == done["started_ns"] / 1e9 and end == done["finished_ns"] / 1e9
            events.append("host_finish")
            if fault == "host_finish":
                raise RuntimeError("host finish failed")
            return {"controlled_workload_verified": False}
        def close(self):
            events.append("host_close")
    monkeypatch.setattr(module, "PeriodicHostObserver", Observer)
    monkeypatch.setattr(module.periodic, "read_job_memory", lambda *a: {})
    monkeypatch.setattr(module.periodic, "step_memory", lambda *a: {})
    monkeypatch.setattr(module.periodic, "interval_point", lambda *a: {})
    monkeypatch.setattr(module.periodic, "evaluate_lineage", lambda p, *a: {"points": len(p)})
    monkeypatch.setattr(module.periodic, "evaluate", lambda p, *a: {"points": len(p)})
    def read(*args):
        number = len(points)
        if fault == "initial" and not number or fault == "final" and number:
            raise ValueError("boundary point failed")
        if fault == "stale" and not number:
            clock[0] += 2.
        t = done["finished_ns"] + 1 if number else 0
        value = dict(thread_affinity=dict(status="observed_within_affinity", errors=[],
            violating_tids=[], initial_tids=[1], final_tids=[1, 2] if number and fault == "descendant" else [1],
            threads=[dict(tid=1, start_ticks=10)], scope="/slurm/job_42/step_0/user",
            started_ns=t, finished_ns=t))
        points.append(value)
        return value
    monkeypatch.setattr(module.periodic, "read_point", read)
    def sleep(duration):
        assert json.loads((directory / "go.json").read_text()) == {"go": True}
        clock[0] += duration
        if fault == "deadline":
            clock[0] = 85831.
        elif clock[0] >= cycles:
            module.save(directory / "done.json", done)
    monkeypatch.setattr(module.time, "sleep", sleep)
    def guard(path):
        assert path == directory and not (directory / "go.json").exists() and not points
        events.append("guard")
        if fault == "guard":
            raise ValueError("release denied")
    return directory, points, events, guard


@pytest.mark.parametrize("cycles", [1, 65, 900])
@pytest.mark.parametrize("code,timed_out", [(0, False), (7, False), (124, True)])
def test_only_two_native_points_and_unchanged_worker(tmp_path, monkeypatch, cycles, code, timed_out):
    directory, points, events, guard = setup(tmp_path, monkeypatch, cycles=cycles, code=code, timed_out=timed_out)
    result = module.measure(["/native"], directory, 42, 32, 128 * 1024**3, 85800, 1., release_guard=guard)
    assert len(points) == 2
    assert result["schema"] == "threadripper_boundary_control_v1"
    assert result["native"]["exit_code"] == code and result["native"]["timed_out"] is timed_out
    assert result["policy"]["native_points"] == 2
    assert result["policy"]["periodic_native_sampling"] is False
    assert result["screening"] == {"points": 2} and result["root_context"] == {"points": 2}
    assert result["scientific_timings_admitted"] is False and result["native_outputs_validated"] is False
    assert events == ["host_initial", "guard", "host_start", "host_finish", "cleanup", "host_close"]
    assert [p.name for p in sorted(directory.glob("point_*.json"))] == ["point_000000.json", "point_000001.json"]
    assert json.loads((directory / "report_finalization.json").read_text())["status"] == "reporting_completed"
    assert not (directory / "lineage_report.json").exists()


@pytest.mark.parametrize("fault", ["guard", "stale", "initial", "final", "descendant", "host_finish",
    "early_exit", "wrapper", "deadline"])
def test_failures_are_retained_without_native_release_or_admission(tmp_path, monkeypatch, fault):
    directory, points, events, guard = setup(tmp_path, monkeypatch, fault=fault)
    with pytest.raises((ValueError, RuntimeError, TimeoutError)):
        module.measure(["/native"], directory, 42, 32, 128 * 1024**3, 85800, 1., release_guard=guard)
    assert (directory / "boundary_failure.json").exists()
    assert not (directory / "boundary_report.json").exists()
    assert json.loads((directory / "release.json").read_text()) == {"release": True}
    if fault in {"guard", "stale", "initial"}:
        assert json.loads((directory / "go.json").read_text()) == {"abort": True}
        assert not (directory / "done.json").exists()
    if fault == "descendant":
        assert json.loads((directory / "native_completion.json").read_text())["errors"]


@pytest.mark.parametrize("key,value", [("job", True), ("cpus", 20), ("memory", 96 * 1024**3),
    ("timeout", 900), ("interval", True), ("monitor", False), ("host_interval", 1), ("guard", "yes")])
def test_wrong_settings_do_not_create_output(tmp_path, monkeypatch, key, value):
    directory, _, _, guard = setup(tmp_path, monkeypatch)
    values = dict(job=42, cpus=32, memory=128 * 1024**3, timeout=85800,
                  interval=1., monitor=True, host_interval=30., guard=guard)
    values[key] = value
    with pytest.raises(ValueError):
        module.measure(["/native"], directory, values["job"], values["cpus"], values["memory"],
            values["timeout"], values["interval"], monitor_host=values["monitor"],
            host_interval_s=values["host_interval"], release_guard=values["guard"])
    assert not directory.exists()

from pathlib import Path
import json
from types import SimpleNamespace

import pytest

from benchmark_tools import measure_threadripper_scaling as collector


@pytest.mark.parametrize("cpus,timeout,cadence", [(32, 85800, 1.), (32, 85800, 1)])
def test_frozen_settings(cpus, timeout, cadence):
    collector.validate(["/bin/sleep", "1"], cpus, timeout, cadence)


@pytest.mark.parametrize("command,cpus,timeout,cadence", [
    (["sleep", "1"], 32, 85800, 1.), ([], 32, 85800, 1.),
    (["/bin/sleep"], 20, 85800, 1.), (["/bin/sleep"], 64, 85800, 1.),
    (["/bin/sleep"], 32, 900, 1.), (["/bin/sleep"], 32, 85800, True),
])
def test_wrong_settings(command, cpus, timeout, cadence):
    with pytest.raises(ValueError):
        collector.validate(command, cpus, timeout, cadence)


def test_step_user_scope_and_retained_gap(monkeypatch):
    monkeypatch.setattr(collector, "read_root_point", lambda *args: {"root": "unchanged"})
    def observe(scope, cpus):
        assert scope == Path("/sys/fs/cgroup/slurm/job_2/step_0/user")
        assert list(cpus) == list(range(32))
        return {"status": "incomplete", "errors": ["thread exited"]}
    monkeypatch.setattr(collector, "observe", observe)
    point = collector.read_point(11, "0::/slurm/job_2/step_0/user/task_0\n", 2, None)
    assert point["root"] == "unchanged"
    assert point["thread_affinity"]["status"] == "incomplete"


def test_missing_user_scope_rejected(monkeypatch):
    monkeypatch.setattr(collector, "read_root_point", lambda *args: {})
    with pytest.raises(ValueError, match="user subtree"):
        collector.read_point(11, "0::/slurm/job_2/step_0/task_0\n", 2, None)


@pytest.mark.parametrize("denied", [None, "budget", "stale", "descendant"])
def test_release_gate_precedes_native_and_denial_aborts(tmp_path, monkeypatch, denied):
    directory = tmp_path / "measurement"
    events = []
    for key, value in {"SLURM_JOB_ID": "42", "SLURM_CPUS_PER_TASK": "64",
                       "SLURM_MEM_PER_NODE": "131072"}.items():
        monkeypatch.setenv(key, value)
    monkeypatch.setattr(collector.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    class Process:
        finished = False
        def poll(self):
            return 0 if self.finished else None
        def wait(self, **kwargs):
            self.finished = True
            events.append("cleanup")
            return 0
    monkeypatch.setattr(collector.subprocess, "Popen", lambda *a, **k: Process())
    monkeypatch.setattr(collector, "wait_file", lambda *a: dict(pid=1,
        cgroup="0::/slurm/job_42/step_0/user/task_0\n", placement={}))
    def host_monitor(*args, sample_fn):
        assert sample_fn is collector.enriched_snapshot
        return SimpleNamespace(observe=lambda: None, summary=lambda *a: {})
    monkeypatch.setattr(collector, "HostMonitor", host_monitor)
    monkeypatch.setattr(collector, "PeriodicHostObserver", lambda host, period:
                        SimpleNamespace(start=lambda: None, close=lambda: None,
                                        finish=lambda *args: host.summary(*args)))
    monkeypatch.setattr(collector, "read_job_memory", lambda *a: {})
    clock = collector.time.monotonic
    offset = [0.]
    monkeypatch.setattr(collector.time, "monotonic", lambda: clock() + offset[0])
    point_count = [0]
    def read_point(*args):
        if denied == "stale":
            offset[0] += 2.
        start = 0 if point_count[0] == 0 else 3
        point_count[0] += 1
        return {"fixture": True, "thread_affinity": dict(status="observed_within_affinity",
            errors=[], violating_tids=[], initial_tids=[1],
            final_tids=[1, 2] if denied == "descendant" and start else [1],
            threads=[dict(tid=1, start_ticks=10)], scope="/slurm/job_42/step_0/user",
            started_ns=start, finished_ns=start)}
    monkeypatch.setattr(collector, "read_point", read_point)
    monkeypatch.setattr(collector, "interval_point", lambda *a: {})
    monkeypatch.setattr(collector, "step_memory", lambda *a: {})
    monkeypatch.setattr(collector, "evaluate_lineage", lambda *a: {})
    monkeypatch.setattr(collector, "evaluate", lambda *a: {})
    def sleep(_):
        assert json.loads((directory / "go.json").read_text()) == {"go": True}
        events.append("native")
        collector.save(directory / "done.json", dict(exit_code=0, started_ns=1, finished_ns=2))
    monkeypatch.setattr(collector.time, "sleep", sleep)
    def guard(path):
        assert path == directory and not (path / "go.json").exists()
        assert not (path / "point_000000.json").exists()
        events.append("guard")
        if denied == "budget":
            raise ValueError("insufficient time")
    def invoke():
        return collector.measure(["/native"], directory, 42, 32, 128*1024**3,
                                 85800, 1., release_guard=guard)
    if denied in {"budget", "stale"}:
        with pytest.raises(ValueError, match="insufficient|stale"):
            invoke()
        assert json.loads((directory / "go.json").read_text()) == {"abort": True}
        assert not (directory / "done.json").exists()
        assert events == ["guard", "cleanup"]
        if denied == "stale":
            assert (directory / "release_freshness_failed.json").exists()
    elif denied == "descendant":
        with pytest.raises(ValueError, match="completion unverified"):
            invoke()
        receipt = json.loads((directory / "native_completion.json").read_text())
        assert receipt["errors"] and receipt["status"] == "native_completion_unverified"
        assert not (directory / "lineage_report.json").exists()
        assert json.loads((directory / "release.json").read_text()) == {"release": True}
        assert events == ["guard", "native", "cleanup"]
    else:
        result = invoke()
        assert result["schema"] == "threadripper_scaling_v5"
        assert result["native_completion"]["status"] == "anchor_only_at_boundaries"
        receipt = json.loads((directory / "report_finalization.json").read_text())
        assert receipt["status"] == "reporting_completed"
        assert receipt["job_id"] == 42
        assert events == ["guard", "native", "cleanup"]

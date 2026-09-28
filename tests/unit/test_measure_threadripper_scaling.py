from pathlib import Path

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

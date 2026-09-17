from copy import deepcopy
from contextlib import contextmanager
from types import SimpleNamespace

import pytest

from benchmark_tools.observe_host_competition import analyze, observe
from benchmark_tools import observe_host_competition as observer


def process(pid, cpu=0., when=1., created=1., group="/other"):
    return {"pid": pid, "created": created, "name": "test", "cgroup": group,
            "user_s": cpu, "system_s": 0., "observed_monotonic_s": when}


def sample(*rows, errors=None):
    return {"processes": list(rows), "errors": errors or []}


def test_own_subtree_and_observer_excluded_but_similar_prefix_not_excluded():
    first = sample(process(1, group="/job/task"), process(2, group="/job/task/child"), process(3),
                   process(4, group="/job/tasks_other"))
    second = deepcopy(first)
    for row in second["processes"]:
        row.update(user_s=2., observed_monotonic_s=3.)
    result = analyze(first, second, "/job/task", 3)
    assert result["status"] == "competing_cpu_observed"
    assert result["sum_observed_foreign_average_cores"] == 1.
    assert [r["pid"] for r in result["persistent_foreign_processes"]] == [4]


def test_small_processes_aggregate_to_competition():
    before = sample(process(1), process(2))
    after = sample(process(1, .2, 2.), process(2, .2, 2.))
    assert analyze(before, after, "/own", 99)["status"] == "competing_cpu_observed"


@pytest.mark.parametrize("change", ["birth", "death", "pid_reuse", "cgroup", "counter", "time", "permission"])
def test_incomplete_observation_cannot_report_clean_window(change):
    first, second = sample(process(1)), sample(process(1, when=2.))
    if change == "birth":
        second["processes"].append(process(2, when=2.))
    elif change == "death":
        second["processes"].clear()
    elif change == "pid_reuse":
        second["processes"][0]["created"] = 2.
    elif change == "cgroup":
        second["processes"][0]["cgroup"] = "/elsewhere"
    elif change == "counter":
        second["processes"][0]["system_s"] = -1.
    elif change == "time":
        second["processes"][0]["observed_monotonic_s"] = 1.
    else:
        second["errors"].append({"pid": 2, "type": "AccessDenied"})
    result = analyze(first, second, "/own", 99)
    assert result["status"] == "inconclusive"
    assert result["controlled_workload_verified"] is False


def test_quiet_sample_is_not_exclusivity_claim():
    result = analyze(sample(process(1)), sample(process(1, when=2.)), "/own", 99)
    assert result["status"] == "no_large_persistent_competitor_observed"
    assert result["controlled_workload_verified"] is False


def test_duplicate_pid_rejected():
    with pytest.raises(ValueError, match="Duplicate"):
        analyze(sample(process(1), process(1, created=2.)), sample(), "/own", 99)


def test_existing_observation_not_overwritten(tmp_path):
    with pytest.raises(FileExistsError):
        observe(1, 1, 1., tmp_path)


@pytest.mark.parametrize("problem", [None, "pid_reuse", "cgroup", "disappeared", "permission"])
def test_snapshot_fresh_identity_check_outside_oneshot(monkeypatch, problem):
    class Process:
        cached = False
        creations = 0

        @contextmanager
        def oneshot(self):
            self.cached = True
            try:
                yield
            finally:
                self.cached = False

        def create_time(self):
            assert self.cached
            self.creations += 1
            return 10.

        def cpu_times(self):
            assert self.cached
            return SimpleNamespace(user=2., system=3.)

        def name(self):
            if problem == "permission":
                raise observer.psutil.AccessDenied(123)
            return "fixture"

        def is_running(self):
            assert not self.cached
            return problem not in {"pid_reuse", "disappeared"}

    proc = Process()
    groups = iter(["/first", "/changed" if problem == "cgroup" else "/first"])
    monkeypatch.setattr(observer.psutil, "pids", lambda: [123])
    monkeypatch.setattr(observer.psutil, "Process", lambda pid: proc)
    monkeypatch.setattr(observer, "membership", lambda pid: next(groups))
    result = observer.snapshot()
    assert proc.creations == 1
    if problem:
        assert not result["processes"] and len(result["errors"]) == 1
    else:
        assert not result["errors"]
        row = result["processes"][0]
        assert (row["pid"], row["created"], row["cgroup"], row["user_s"], row["system_s"]) == (123, 10., "/first", 2., 3.)

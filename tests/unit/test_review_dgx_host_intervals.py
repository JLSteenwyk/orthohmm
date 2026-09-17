import pytest

from benchmark_tools.review_dgx_host_intervals import unmatched_inventory


def row(processes, unmatched=None):
    return {"snapshot": {"processes": processes, "errors": []},
            "interval": None if unmatched is None else {
                "status": "inconclusive", "uncertain_processes": [], "unmatched_foreign_processes": unmatched}}


def test_pid_reuse_preserves_distinct_identity_directions():
    before = dict(pid=5, created=1, name="kworker/1", cgroup="/")
    after = dict(pid=5, created=2, name="other", cgroup="/user")
    data = unmatched_inventory([row([before]), row([after], [dict(pid=5, created=1), dict(pid=5, created=2)])])
    assert data["unmatched_identity_events"] == 2
    assert data["kworker_named_identity_events"] == 1
    assert data["unmatched_events"] == [
        dict(name="kworker/1", cgroup="/", direction="disappeared", count=1),
        dict(name="other", cgroup="/user", direction="appeared", count=1)]


def test_matched_identity_cannot_be_labeled_unmatched():
    process = dict(pid=5, created=1, name="worker", cgroup="/")
    with pytest.raises(ValueError, match="exclusive"):
        unmatched_inventory([row([process]), row([process], [dict(pid=5, created=1)])])


def test_observation_errors_are_not_silently_dropped():
    with pytest.raises(ValueError, match="separate review"):
        unmatched_inventory([{"observation_error": "PermissionError"}])


def test_duplicate_identity_rejected():
    process = dict(pid=5, created=1, name="worker", cgroup="/")
    with pytest.raises(ValueError, match="Duplicate"):
        unmatched_inventory([row([process, process])])

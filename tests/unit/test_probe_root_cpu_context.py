import copy

import pytest

from benchmark_tools import probe_root_cpu_context as probe


@pytest.fixture
def tree(tmp_path):
    root, proc = tmp_path / "cgroup", tmp_path / "proc"
    for name in probe.SCOPES:
        path = root / name.lstrip("/")
        path.mkdir(parents=True, exist_ok=True)
        (path / "cpu.stat").write_text("usage_usec 100\nuser_usec 40\nsystem_usec 60\n")
    (root / "cgroup.procs").write_text("2\n9\n")
    (proc / "sys/kernel/random").mkdir(parents=True)
    (proc / "sys/kernel/random/boot_id").write_text("boot\n")
    (proc / "stat").write_text("cpu  100 0 100 500 0 1 2 0 0 0\n")
    return root, proc


def test_snapshot_and_signed_complement(tree):
    root, proc = tree
    left = probe.snapshot(root, proc)
    for name, value in zip(probe.SCOPES, [110, 120, 105, 101]):
        (root / name.lstrip("/") / "cpu.stat").write_text(f"usage_usec {value}\n")
    (root / "cgroup.procs").write_text("2\n10\n")
    (proc / "stat").write_text("cpu  101 0 101 501 0 1 2 0 0 0\n")
    result = probe.compare(left, probe.snapshot(root, proc))
    assert result["root_minus_system_cpu_usec"] == -10
    assert result["root_minus_three_named_children_cpu_usec"] == -16
    assert result["scope_cpu_usec"]["/user.slice"] == 5
    assert result["enclosing_host_category_ticks"]["user"] == 1
    assert result["observed_root_membership_changed"] is True
    assert result["root_membership_snapshots"] == [[2, 9], [2, 9], [2, 10], [2, 10]]
    assert result["scientific_timings_admitted"] is False
    assert result["environmental_validity_established"] is False


def test_process_churn_does_not_invalidate_stable_scope_counters(tree, monkeypatch):
    root, proc = tree
    original = probe.read_counter

    def read_counter(path, scope):
        result = original(path, scope)
        (root / "cgroup.procs").write_text("5\n")
        return result

    monkeypatch.setattr(probe, "read_counter", read_counter)
    result = probe.snapshot(root, proc)
    assert probe.pids(result["members_before"]["raw"]) == [2, 9]
    assert probe.pids(result["members_after"]["raw"]) == [5]


def test_missing_scope_retains_failure_evidence(tree):
    root, proc = tree
    (root / "user.slice/cpu.stat").unlink()
    with pytest.raises(probe.RootContextError) as error:
        probe.snapshot(root, proc)
    assert error.value.evidence["status"] == "invalid_root_cpu_context"
    assert len(error.value.evidence["rows"]) == 2
    assert error.value.evidence["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["boot", "identity", "alias", "rows", "counter", "ticks", "order", "members", "host"])
def test_malformed_context_rejected(tree, fault):
    point = probe.snapshot(*tree)
    if fault == "boot":
        point["boot_after"] = "other"
    elif fault == "identity":
        point["identities_after"]["/"][1] += 1
    elif fault == "alias":
        for key in ("identities_before", "identities_after"):
            point[key]["/init.scope"] = point[key]["/"]
    elif fault == "rows":
        point["rows"].reverse()
    elif fault == "counter":
        point["rows"][0]["raw"] = "usage_usec -1\n"
    elif fault == "ticks":
        point["ticks"] = True
    elif fault == "order":
        point["rows"][0]["started_ns"] = 0
    elif fault == "members":
        point["members_after"]["raw"] = "-2\n"
    else:
        point["host_after"]["raw"] = "cpu  1 0 1 1 0 0 0 0 0 0\n"
    with pytest.raises(ValueError):
        probe.validate(point)


@pytest.mark.parametrize("fault", ["boot", "identity", "overlap", "decrease"])
def test_invalid_comparison_rejected(tree, fault):
    left, right = probe.snapshot(*tree), probe.snapshot(*tree)
    if fault == "boot":
        right["boot_before"] = right["boot_after"] = "other"
    elif fault == "identity":
        for key in ("identities_before", "identities_after"):
            right[key]["/"][1] += 100
    elif fault == "overlap":
        right = copy.deepcopy(left)
    else:
        right["rows"][0]["raw"] = "usage_usec 1\n"
    with pytest.raises(ValueError):
        probe.compare(left, right)


def test_membership_normalizes_duplicates_without_claiming_continuous_identity():
    assert probe.pids("9\n2\n9\n") == [2, 9]
    assert probe.pids("") == []

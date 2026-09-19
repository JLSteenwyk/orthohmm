import copy
import pytest

from benchmark_tools import probe_cgroup_lineage as probe


TARGET = "/system.slice/slurm.scope/job_1"


def point(offset=0, counters=(100, 80, 50, 30)):
    names = probe.scopes(TARGET)
    ids = {name: [1, i + 10] for i, name in enumerate(names)}
    return dict(status="aggregate_lineage_snapshot", target=TARGET,
        boot_before="boot", boot_after="boot", identities_before=ids,
        identities_after=copy.deepcopy(ids), rows=[dict(scope=name,
            started_ns=offset + i * 10, finished_ns=offset + i * 10 + 1,
            raw=f"usage_usec {value}\nuser_usec 0\nsystem_usec 0\n")
            for i, (name, value) in enumerate(zip(names, counters))])


def test_aggregate_complements_not_overlapping_sum():
    result = probe.compare(point(), point(100, (200, 160, 110, 80)))
    assert [r["cpu_usec"] for r in result["aggregate_deltas"]] == [100, 80, 60, 50]
    assert [r["signed_cpu_usec"] for r in result["signed_complements"]] == [20, 20, 10]
    assert result["root_minus_target_cpu_usec"] == 50
    assert not result["scientific_timings_admitted"]
    assert not result["controlled_workload_verified"]


def test_negative_complement_retained():
    result = probe.compare(point(), point(100, (120, 110, 90, 80)))
    assert result["root_minus_target_cpu_usec"] == -30
    assert [r["signed_cpu_usec"] for r in result["signed_complements"]] == [-10, -10, -10]


@pytest.mark.parametrize("target", ["/", "relative", "/a/../b", "/a//b", "/a/"])
def test_invalid_target(target):
    with pytest.raises(ValueError):
        probe.scopes(target)


@pytest.mark.parametrize("change", ["status", "boot", "identity", "missing", "alias", "rows",
                                    "duplicate", "order", "bool_time", "counter"])
def test_invalid_snapshot(change):
    value = point()
    if change == "status":
        value["status"] = "invalid_lineage_snapshot"
    elif change == "boot":
        value["boot_after"] = "other"
    elif change == "identity":
        value["identities_after"][TARGET][1] += 1
    elif change == "missing":
        for key in ("identities_before", "identities_after"):
            value[key].pop("/")
    elif change == "alias":
        for key in ("identities_before", "identities_after"):
            value[key][TARGET] = value[key]["/"]
    elif change == "rows":
        value["rows"].pop()
    elif change == "duplicate":
        value["rows"][1] = value["rows"][0]
    elif change == "order":
        value["rows"][1]["started_ns"] = 0
    elif change == "bool_time":
        value["rows"][0]["started_ns"] = False
    else:
        value["rows"][0]["raw"] = "usage_usec -1\n"
    with pytest.raises(ValueError):
        probe.validate(value)


@pytest.mark.parametrize("change", ["boot", "identity", "overlap", "decrease"])
def test_invalid_interval(change):
    left, right = point(), point(100, (200, 160, 110, 80))
    if change == "boot":
        right["boot_before"] = right["boot_after"] = "other"
    elif change == "identity":
        for key in ("identities_before", "identities_after"):
            right[key][TARGET][1] += 1
    elif change == "overlap":
        right = point(10, (200, 160, 110, 80))
    else:
        right["rows"][-1]["raw"] = "usage_usec 0\n"
    with pytest.raises(ValueError):
        probe.compare(left, right)


def hierarchy(tmp_path):
    root = tmp_path / "cgroups"
    for name in probe.scopes(TARGET):
        directory = root / name.lstrip("/")
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "cpu.stat").write_text("usage_usec 100\n")
    boot = tmp_path / "boot"
    boot.write_text("boot\n")
    return root, boot


def test_control_callback_runs_after_each_retained_counter(tmp_path):
    root, boot = hierarchy(tmp_path)
    visited = []
    result = probe.snapshot(root, TARGET, boot, after_read=visited.append)
    assert visited == [r["scope"] for r in result["rows"]] == probe.scopes(TARGET)


def test_failed_control_callback_retains_partial_snapshot(tmp_path):
    root, boot = hierarchy(tmp_path)

    def fail(scope):
        raise ValueError("injected failure")

    with pytest.raises(probe.LineageSnapshotError) as error:
        probe.snapshot(root, TARGET, boot, after_read=fail)
    assert len(error.value.evidence["rows"]) == 1
    assert error.value.evidence["rows"][0]["scope"] == "/"


def test_sibling_churn_during_read_does_not_invalidate_lineage(tmp_path, monkeypatch):
    root, boot = hierarchy(tmp_path)
    original = probe.read_counter
    def read(directory, scope):
        sibling = root / "system.slice/transient.service"
        if scope == "/":
            sibling.mkdir()
        if scope == TARGET:
            sibling.rmdir()
        return original(directory, scope)
    monkeypatch.setattr(probe, "read_counter", read)
    result = probe.snapshot(root, TARGET, boot)
    assert len(result["rows"]) == 4
    assert not (root / "system.slice/transient.service").exists()
    probe.validate(result)


def test_lineage_change_retains_failed_evidence(tmp_path, monkeypatch):
    root, boot = hierarchy(tmp_path)
    original = probe.read_counter
    def read(directory, scope):
        row = original(directory, scope)
        if scope == TARGET:
            target = root / TARGET.lstrip("/")
            target.rename(target.with_name("old_job"))
            target.mkdir()
        return row
    monkeypatch.setattr(probe, "read_counter", read)
    with pytest.raises(probe.LineageSnapshotError) as caught:
        probe.snapshot(root, TARGET, boot)
    assert len(caught.value.evidence["rows"]) == 4
    assert caught.value.evidence["status"] == "invalid_lineage_snapshot"


def test_symlink_rejected(tmp_path):
    root, boot = hierarchy(tmp_path)
    target = root / TARGET.lstrip("/")
    target.rename(target.with_name("old_job"))
    target.symlink_to(target.with_name("old_job"), target_is_directory=True)
    with pytest.raises(probe.LineageSnapshotError):
        probe.snapshot(root, TARGET, boot)


def test_unreadable_counter_preserves_partial_rows(tmp_path):
    root, boot = hierarchy(tmp_path)
    (root / TARGET.lstrip("/") / "cpu.stat").unlink()
    with pytest.raises(probe.LineageSnapshotError) as caught:
        probe.snapshot(root, TARGET, boot)
    assert len(caught.value.evidence["rows"]) == 3
    assert caught.value.evidence["error_type"] == "FileNotFoundError"
    assert not caught.value.evidence["scientific_timings_admitted"]


def test_boot_change_during_snapshot_rejected(tmp_path, monkeypatch):
    root, boot = hierarchy(tmp_path)
    original = probe.read_counter
    def read(directory, scope):
        boot.write_text("other\n")
        return original(directory, scope)
    monkeypatch.setattr(probe, "read_counter", read)
    with pytest.raises(probe.LineageSnapshotError, match="boot identity"):
        probe.snapshot(root, TARGET, boot)

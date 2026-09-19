from copy import deepcopy

import pytest

from benchmark_tools import probe_native_pressure as module
from tests.unit.test_probe_dgx_cpu_hierarchy import point as hierarchy_point


def raw(some, full):
    return (f"some avg10=0.00 avg60=0.00 avg300=0.00 total={some}\n"
            f"full avg10=0.00 avg60=0.00 avg300=0.00 total={full}\n")


def point(index):
    old = hierarchy_point(index)
    pressure = {r: dict(raw=raw(index*100, index*40), totals=dict(some=index*100, full=index*40),
                       started_ns=index*1_000_000_000+20+i*10,
                       finished_ns=index*1_000_000_000+21+i*10)
                for i, r in enumerate(module.RESOURCES)}
    return dict(host=old["host"], native_membership=old["native_membership"], ticks=100,
                scope="/job_123/step_0", scope_identity=[32, 100], pressure=pressure)


def test_native_cpu_full_is_meaningful_not_system_cpu_full():
    result = module.compare(point(0), point(1), 123)
    assert result["native_stall_usec"]["cpu"] == dict(some=100, full=40)
    assert result["scientific_timings_admitted"] is False
    assert result["controlled_workload_verified"] is False


@pytest.mark.parametrize("fault", ["inode", "boot", "scope", "order", "raw", "missing", "counter", "overlap", "batch"])
def test_invalid_native_evidence_rejected(fault):
    left, right = point(0), point(1)
    if fault == "inode":
        right["scope_identity"][1] += 1
    elif fault == "boot":
        for h in right["host"]:
            h["raw"]["boot_id"] = "different"
    elif fault == "scope":
        right["scope"] += "/task_0"
    elif fault == "order":
        right["pressure"]["cpu"]["started_ns"] = 0
    elif fault == "raw":
        right["pressure"]["cpu"]["totals"]["some"] = 200
    elif fault == "missing":
        del right["pressure"]["memory"]
    elif fault == "counter":
        left["pressure"]["cpu"].update(raw=raw(200, 100), totals=dict(some=200, full=100))
    elif fault == "overlap":
        right = deepcopy(left)
    else:
        right["native_membership"] = "0::/job_123/step_batch/task_0\n"
    with pytest.raises(ValueError):
        module.compare(left, right, 123)


def test_reads_pressure_at_step_not_task_scope(tmp_path, monkeypatch):
    p = point(0)
    proc = tmp_path / "proc"
    (proc / "42").mkdir(parents=True)
    (proc / "42/cgroup").write_text(p["native_membership"])
    groups = tmp_path / "groups"
    directory = groups / "job_123/step_0"
    directory.mkdir(parents=True)
    for resource in module.RESOURCES:
        (directory / f"{resource}.pressure").write_text(raw(0, 0))
    times = iter([20, 21, 30, 31, 40, 41])
    monkeypatch.setattr(module.time, "monotonic_ns", lambda: next(times))
    monkeypatch.setattr(module.os, "sysconf", lambda _: 100)
    hosts = iter(p["host"])
    result = module.read_point(42, p["native_membership"], 123, proc_root=proc, group_root=groups,
                               host_reader=lambda: next(hosts))
    assert result["pressure"] == p["pressure"]
    assert result["scope"] == p["scope"]
    assert result["scope_identity"] == [directory.stat().st_dev, directory.stat().st_ino]


def test_unknown_or_changed_process_rejected(tmp_path):
    (tmp_path / "42").mkdir()
    (tmp_path / "42/cgroup").write_text("changed")
    with pytest.raises(ValueError, match="membership changed"):
        module.read_point(42, "expected", 123, proc_root=tmp_path)

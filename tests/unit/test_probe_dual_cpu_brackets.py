from copy import deepcopy

import pytest

from benchmark_tools import probe_dual_cpu_brackets as module
from tests.unit.test_measure_native_hierarchy_step import evidence
from tests.unit.test_measure_native_frontier_step import extend
from tests.unit.test_frontier_pressure_integration import add_pressure


def paired(evidence):
    points, done = evidence
    extend(points, done)
    add_pressure(points)
    return points


def test_same_native_scope_and_read_times(evidence):
    points = paired(evidence)
    before = deepcopy(points)
    result = module.compare(*points, 21816)
    assert result["outer"]["native_cpu_s"] == result["narrow"]["native_cpu_s"]
    assert result["outer"]["wall_s"] == result["narrow"]["wall_s"]
    assert result["outer"]["outer_read_overhang_s"] > result["narrow"]["outer_read_overhang_s"]
    assert result["native_pressure"]["native_stall_usec"]["cpu"]["some"] == 100
    assert result["frontier"]["outside_target_frontier_cpu_s"] == pytest.approx(.02)
    assert points == before
    assert result["scientific_timings_admitted"] is False
    assert result["controlled_workload_verified"] is False


def test_reader_requests_pressure_and_preserves_failure_path(evidence, monkeypatch, tmp_path):
    point = paired(evidence)[0]
    failure = tmp_path / "failed.json"
    def read(pid, membership, job, path, *, native_pressure):
        assert (pid, membership, job, path, native_pressure) == (42, "membership", 21816, failure, True)
        return point
    monkeypatch.setattr(module, "read_frontier_point", read)
    assert module.read_point(42, "membership", 21816, failure) is point


@pytest.mark.parametrize("fault", ["missing_pressure", "pressure_identity", "frontier_identity", "boot", "counter", "time"])
def test_invalid_paired_evidence_rejected(evidence, fault):
    left, right = paired(evidence)
    if fault == "missing_pressure":
        del right["native_pressure"]
    elif fault == "pressure_identity":
        right["native_pressure"]["scope_identity"][1] += 1
    elif fault == "frontier_identity":
        for inv in ("inventory_before", "inventory_after"):
            right["frontier"][inv]["identities"]["/user.slice"][1] += 1
    elif fault == "boot":
        right["hierarchy_host_after"]["raw"]["boot_id"] = "other"
    elif fault == "counter":
        right["hierarchy_host_after"]["cpu_ticks"]["user"] += 1
    else:
        right["hierarchy_host_after"]["finished_monotonic_ns"] += 1000000000
    with pytest.raises(ValueError):
        module.compare(left, right, 21816)


def test_whole_window_opt_out_does_not_bypass_identity(evidence):
    left, right = paired(evidence)
    expected = module.compare(left, right, 21816)
    assert module.compare(left, right, 21816, enforce_gap=False) == expected
    right["native_pressure"]["scope_identity"][1] += 1
    with pytest.raises(ValueError):
        module.compare(left, right, 21816, enforce_gap=False)

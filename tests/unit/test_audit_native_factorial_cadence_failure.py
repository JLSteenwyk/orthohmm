from copy import deepcopy

import pytest

from benchmark_tools.audit_native_factorial_cadence_failure import ERROR, cadence
from tests.unit.test_probe_interval_cpu import point


def retime(value, ns):
    value = deepcopy(value)
    delta = ns - value["host"][0]["started_monotonic_ns"]
    for host in value["host"]:
        for key in ("started_monotonic_ns", "finished_monotonic_ns"):
            host[key] += delta
    value["native_read_ns"] = [n + delta for n in value["native_read_ns"]]
    return value


def test_valid_cadence_never_admits_timing_outputs_or_successor():
    report = cadence((point(i) for i in range(4)), 123)
    assert report["points"] == 4 and report["intervals"] == 3
    assert report["failures"] == [] and report["failure_counts"] == {}
    for key in ("scientific_timings_admitted", "native_outputs_validated",
                "next_identity_authorized", "full_resource_replay"):
        assert report[key] is False


@pytest.mark.parametrize("seconds,failed", [(.5, False), (1.5, False), (.499, True), (1.501, True)])
def test_exact_unchanged_cadence_bounds(seconds, failed):
    report = cadence([point(0), retime(point(1), int(seconds * 1e9))], 123)
    assert report["unchanged_cadence_bounds_s"] == [.5, 1.5]
    assert bool(report["failures"]) is failed
    if failed:
        assert report["failures"][0]["error"] == ERROR
        assert report["failures"][0]["wall_s"] == seconds


def test_all_slow_and_catchup_intervals_preserved_without_interpolation():
    report = cadence([point(0), retime(point(1), 1_700_000_000), point(2), point(3)], 123)
    assert report["points"] == 4
    assert report["failure_counts"] == {ERROR: 2}
    assert [(f["left_index"], f["right_index"], f["wall_s"]) for f in report["failures"]] == [
        (0, 1, 1.7), (1, 2, .3)]


def test_non_cadence_pair_rejection_is_not_hidden():
    right = point(1)
    right["native_cpu_stat"] = "usage_usec 0\n"
    left = point(0)
    left["native_cpu_stat"] = "usage_usec 1000000\n"
    report = cadence([left, right], 123)
    assert report["failure_counts"] == {"Native CPU counter decreased": 1}


def test_invalid_point_still_rejected():
    right = point(1)
    right["native_cpu_scope"] = "/job_999/step_0"
    with pytest.raises(ValueError):
        cadence([point(0), right], 123)


@pytest.mark.parametrize("count", [0, 1])
def test_incomplete_inventory_rejected(count):
    with pytest.raises(ValueError, match="at least two"):
        cadence([point(i) for i in range(count)], 123)

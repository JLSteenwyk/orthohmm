from copy import deepcopy

import pytest

from benchmark_tools.diagnose_pressure_brackets import compare
from benchmark_tools.measure_native_hierarchy_step import evaluate
from benchmark_tools.probe_host_counters import parse_cpu
from tests.unit.test_measure_native_hierarchy_step import evidence
from tests.unit.test_measure_native_frontier_step import extend


def build(evidence):
    points, done = evidence
    extend(points, done)
    return dict(points=points, job_id=21816, screening=evaluate(points, done, 21816))


def test_replay_keeps_same_native_counters_and_does_not_mutate(evidence):
    report = build(evidence)
    before = deepcopy(report)
    result = compare(report)
    assert report == before
    assert result["intervals"] == len(report["points"]) - 1
    assert result["added_host_cpu_s"]["maximum"] == 0
    assert result["added_read_overhang_s"]["minimum"] == pytest.approx(.001)
    assert sum(result["transitions"].values()) == result["intervals"]
    assert result["scientific_timings_admitted"] is False


def test_extra_host_ticks_are_isolated_without_native_change(evidence):
    points, done = evidence
    extend(points, done)
    for point in points:
        host = point["host"][1]
        lines = host["raw"]["proc_stat"].splitlines()
        fields = lines[0].split()
        fields[1] = str(int(fields[1]) + 30)
        lines[0] = " ".join(fields)
        host["raw"]["proc_stat"] = "\n".join(lines) + "\n"
        host["cpu_ticks"] = parse_cpu(host["raw"]["proc_stat"])
    report = dict(points=points, job_id=21816, screening=evaluate(points, done, 21816))
    result = compare(report)
    assert result["added_host_cpu_s"]["minimum"] == pytest.approx(.3)
    assert result["added_host_cpu_s"]["maximum"] == pytest.approx(.3)


@pytest.mark.parametrize("fault", ["time", "boot", "counter", "stored", "flags", "missing"])
def test_invalid_evidence_rejected(evidence, fault):
    report = build(evidence)
    earlier = report["points"][0]["hierarchy_host_after"]
    if fault == "time":
        earlier["finished_monotonic_ns"] += 10000000000
    elif fault == "boot":
        earlier["raw"]["boot_id"] = "other"
    elif fault == "counter":
        earlier["cpu_ticks"]["user"] += 1
    elif fault == "stored":
        report["screening"]["original_threshold_screen"]["intervals"][0]["native_cpu_s"] += 1
    elif fault == "flags":
        report["screening"]["original_threshold_screen"]["flagged_intervals"] = [-1]
    else:
        report["points"] = []
    with pytest.raises(ValueError):
        compare(report)

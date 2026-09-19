from copy import deepcopy
import inspect
import json

import pytest

from benchmark_tools import measure_native_lineage_step as module
from benchmark_tools import measure_native_dual_bracket_step as original
from benchmark_tools.probe_cgroup_lineage import scopes
from tests.unit.test_measure_native_hierarchy_step import evidence, setup_measure
from tests.unit.test_measure_native_frontier_step import extend
from tests.unit.test_frontier_pressure_integration import add_pressure


def extend_lineage(evidence):
    points, done = evidence
    extend(points, done)
    add_pressure(points)
    for index, point in enumerate(points):
        frontier = point.pop("frontier")
        names = scopes(frontier["target"])
        start = frontier["root"][0]["started_ns"]
        ids = {name: [1, i + 100] for i, name in enumerate(names)}
        point["schema"] = "native_lineage_v1"
        point["lineage"] = dict(status="aggregate_lineage_snapshot", target=frontier["target"],
            boot_before=frontier["boot_id"], boot_after=frontier["boot_id"],
            identities_before=deepcopy(ids), identities_after=deepcopy(ids),
            rows=[dict(scope=name, started_ns=start + i * 2, finished_ns=start + i * 2 + 1,
                raw=point["parent"][1]["raw"] if name == frontier["target"]
                    else f"usage_usec {1000000 + index * 100000}\n") for i, name in enumerate(names)])
    return points, done


def test_lifecycle_matches_existing_collector():
    expected = inspect.getsource(original.measure).replace("final dual observation", "final lineage observation")
    expected = expected.replace("dual_bracket_report.json", "lineage_report.json")
    expected = expected.replace("Complete-command dual-bracket engineering measurement", "Complete-command lineage engineering measurement")
    expected = expected.replace("host_interval_s=30.):", "host_interval_s=30., *, point_reader=None):")
    expected = expected.replace("    validate(command, cpus, timeout_s, interval_s)\n",
        "    validate(command, cpus, timeout_s, interval_s)\n"
        "    reader = read_point if point_reader is None else point_reader\n"
        "    if not callable(reader):\n"
        "        raise ValueError(\"Point reader must be callable\")\n")
    expected = expected.replace("points = [read_point(", "points = [reader(")
    expected = expected.replace("points.append(read_point(", "points.append(reader(")
    assert inspect.getsource(module.measure) == expected


def test_same_native_sample_and_original_thresholds(evidence):
    points, done = extend_lineage(evidence)
    before = deepcopy(points)
    result = module.evaluate(points, done, 21816)
    assert result["original_screening"] == module.hierarchy_evaluate(points, done, 21816)
    interval = result["intervals"][0]
    assert interval["outer"]["native_cpu_s"] == interval["narrow"]["native_cpu_s"]
    assert interval["outer"]["wall_s"] == interval["narrow"]["wall_s"]
    assert interval["outer"]["outer_read_overhang_s"] > interval["narrow"]["outer_read_overhang_s"]
    assert points == before
    assert not result["scientific_timings_admitted"]
    assert "frontier" not in interval and "lineage" in interval


@pytest.mark.parametrize("fault", ["schema", "frontier", "target", "boot", "counter", "identity",
    "lineage_early", "lineage_late", "pressure_missing", "pressure_scope", "pressure_identity", "pressure_late", "narrow_time"])
def test_corrupt_observation_rejected(evidence, fault):
    points, done = extend_lineage(evidence)
    point = points[-1]
    if fault == "schema":
        point["schema"] = "old"
    elif fault == "frontier":
        point["frontier"] = {}
    elif fault == "target":
        point["lineage"]["target"] = "/other"
    elif fault == "boot":
        point["lineage"]["boot_before"] = point["lineage"]["boot_after"] = "other"
    elif fault == "counter":
        point["lineage"]["rows"][-1]["raw"] = "usage_usec 0\n"
    elif fault == "identity":
        for key in ("identities_before", "identities_after"):
            point["lineage"][key]["/"][1] += 1
    elif fault == "lineage_early":
        point["lineage"]["rows"][0]["started_ns"] = 0
    elif fault == "lineage_late":
        point["lineage"]["rows"][-1]["finished_ns"] = point["native_pressure"]["host"][0]["started_monotonic_ns"] + 1
    elif fault == "pressure_missing":
        del point["native_pressure"]
    elif fault == "pressure_scope":
        point["native_pressure"]["native_membership"] = "0::/different"
    elif fault == "pressure_identity":
        point["native_pressure"]["scope_identity"][1] += 1
    elif fault == "pressure_late":
        point["native_pressure"]["host"][1]["finished_monotonic_ns"] = point["host"][1]["started_monotonic_ns"] + 1
    else:
        point["hierarchy_host_after"]["finished_monotonic_ns"] += 1000000000
    with pytest.raises((ValueError, KeyError)):
        module.evaluate(points, done, 21816)


@pytest.mark.parametrize("exit_code", [0, 7, 124])
def test_native_status_memory_and_release_preserved(tmp_path, monkeypatch, evidence, exit_code):
    import tests.unit.test_measure_native_hierarchy_step as fixture
    points, done = extend_lineage(evidence)
    monkeypatch.setattr(fixture, "module", module)
    directory = setup_measure(tmp_path, monkeypatch, evidence, exit_code)
    pending = iter(points)
    monkeypatch.setattr(module, "read_point", lambda *args: next(pending))
    result = module.measure(["/usr/bin/true"], directory, 21816, 20, 96 * 1024**3, 60, 1.)
    assert result["native"]["exit_code"] == exit_code
    assert result["native"]["timed_out"] == (exit_code == 124)
    assert result["step_memory"]["scope"] == module.interval_point(points[-1], 21816)["native_cpu_scope"]
    assert json.loads((directory / "lineage_report.json").read_text()) == result
    assert (directory / "release.json").exists()


def test_failure_releases_worker_and_does_not_write_success(tmp_path, monkeypatch, evidence):
    import tests.unit.test_measure_native_hierarchy_step as fixture
    extend_lineage(evidence)
    monkeypatch.setattr(fixture, "module", module)
    directory = setup_measure(tmp_path, monkeypatch, evidence)
    def fail(*args):
        raise ValueError("invalid observation")
    monkeypatch.setattr(module, "read_point", fail)
    with pytest.raises(ValueError):
        module.measure(["/usr/bin/true"], directory, 21816, 20, 96 * 1024**3, 60, 1.)
    assert (directory / "release.json").exists() and (directory / "go.json").exists()
    assert not (directory / "lineage_report.json").exists()


def test_reader_and_migration_failure_retention(tmp_path, monkeypatch, evidence):
    initial = deepcopy(evidence[0][0])
    points, done = extend_lineage(evidence)
    expected = points[0]
    monkeypatch.setattr(module, "read_hierarchy", lambda *args: deepcopy(initial))
    monkeypatch.setattr(module, "lineage_snapshot", lambda *args: deepcopy(expected["lineage"]))
    monkeypatch.setattr(module, "read_pressure", lambda *args: deepcopy(expected["native_pressure"]))
    monkeypatch.setattr(module, "host_snapshot", lambda: deepcopy(expected["host"][1]))
    original_read = module.Path.read_text
    monkeypatch.setattr(module.Path, "read_text", lambda p, *a, **kw:
        initial["native_membership"] if str(p) == "/proc/123/cgroup" else original_read(p, *a, **kw))
    assert module.read_point(123, initial["native_membership"], 21816) == expected
    failure = tmp_path / "failed.json"
    with pytest.raises(ValueError, match="changed scope"):
        module.read_point(123, "other", 21816, failure)
    saved = json.loads(failure.read_text())
    assert saved["status"] == "invalid_native_lineage_observation"
    assert saved["preceding_point"]["lineage"] == expected["lineage"]


def test_partial_lineage_error_retained(tmp_path, monkeypatch, evidence):
    monkeypatch.setattr(module, "read_hierarchy", lambda *args: deepcopy(evidence[0][0]))
    error = module.LineageSnapshotError(ValueError("identity changed"), {"rows": []})
    def fail(*args):
        raise error
    monkeypatch.setattr(module, "lineage_snapshot", fail)
    failure = tmp_path / "failed.json"
    with pytest.raises(module.LineageSnapshotError):
        module.read_point(123, "membership", 21816, failure)
    assert json.loads(failure.read_text())["partial_lineage"] == error.evidence

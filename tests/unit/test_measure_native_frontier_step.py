from copy import deepcopy
import inspect
import json

import pytest

from benchmark_tools import measure_native_frontier_step as module
from benchmark_tools import measure_native_hierarchy_step as original
from tests.unit.test_measure_native_hierarchy_step import evidence, setup_measure


def extend(points, done):
    for index, point in enumerate(points):
        old = point["host"][1]
        point["hierarchy_host_after"] = deepcopy(old)
        point["host"][1] = deepcopy(old)
        for key in ("started_monotonic_ns", "finished_monotonic_ns"):
            point["host"][1][key] += 1_000_000
        start = old["finished_monotonic_ns"] + 1
        target = point["parent"][0]["scope"]
        names = sorted([target, "/user.slice"])
        inventory = dict(identities={name: [1, i] for i, name in enumerate(names)},
                         ancestor_direct_process_counts={"/": 0})
        def row(scope, position, raw):
            return dict(scope=scope, started_ns=start+position*10, finished_ns=start+position*10+1, raw=raw)
        point["frontier"] = dict(target=target, boot_id=old["raw"]["boot_id"].strip(),
            inventory_before=deepcopy(inventory), inventory_after=deepcopy(inventory),
            root=[row("/", 0, f"usage_usec {1000000+index*100000}\n"),
                  row("/", 3, f"usage_usec {1000000+index*100000}\n")],
            rows=[row(name, i+1, point["parent"][1]["raw"] if name == target
                      else f"usage_usec {index*20000}\n") for i, name in enumerate(names)])
    # The worker starts only after the expanded pre-command observation.
    for key in ("started_monotonic_ns", "finished_monotonic_ns"):
        done["snapshots"][0][key] += 2_000_000
    done["started_ns"] += 2_000_000
    return points


def test_collector_body_preserves_worker_lifecycle():
    expected = inspect.getsource(original.measure).replace("read_hierarchy(", "reader(")
    expected = expected.replace('host_interval_s=30.):', 'host_interval_s=30., *, native_pressure=False):\n    reader = partial(read_frontier_point, native_pressure=True) if native_pressure else read_frontier_point')
    expected = expected.replace('ready["cgroup"], job_id)', 'ready["cgroup"], job_id, directory / "failed_point.json")')
    expected = expected.replace("hierarchy_report.json", "frontier_report.json")
    expected = expected.replace("final hierarchy observation", "final frontier observation")
    expected = expected.replace("Complete-command hierarchy engineering test", "Complete-command frontier engineering test")
    assert inspect.getsource(module.measure) == expected


def test_same_original_threshold_screen_with_additional_descriptive_fields(evidence):
    points, done = evidence
    extend(points, done)
    expected = original.evaluate(points, done, 21816)
    result = module.evaluate(points, done, 21816)
    assert all(result[key] == value for key, value in expected.items())
    assert result["frontier_intervals"][0]["outside_target_frontier_cpu_s"] == pytest.approx(.02)
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["target", "boot", "time", "missing", "decrease"])
def test_inconsistent_frontier_rejected(evidence, fault):
    points, done = evidence
    extend(points, done)
    frontier = points[0]["frontier"]
    if fault == "target":
        frontier["target"] = "/user.slice"
    elif fault == "boot":
        frontier["boot_id"] = "changed"
    elif fault == "time":
        frontier["root"][0]["started_ns"] = 0
    elif fault == "missing":
        frontier["rows"].pop()
    else:
        next(row for row in frontier["rows"] if row["scope"] == frontier["target"])["raw"] = "usage_usec 0\n"
    with pytest.raises(ValueError):
        module.evaluate(points, done, 21816)


@pytest.mark.parametrize("exit_code", [0, 7, 124])
def test_native_status_preserved(tmp_path, monkeypatch, evidence, exit_code):
    points, done = evidence
    extend(points, done)
    # Reuse the existing mock worker contract with this collector's globals.
    import tests.unit.test_measure_native_hierarchy_step as tests
    monkeypatch.setattr(tests, "module", module)
    directory = setup_measure(tmp_path, monkeypatch, evidence, exit_code)
    pending = iter(points)
    monkeypatch.setattr(module, "read_frontier_point", lambda *args: next(pending))
    result = module.measure(["/usr/bin/true"], directory, 21816, 20, 96*1024**3, 60, 1.)
    assert result["native"]["exit_code"] == exit_code
    assert result["native"]["timed_out"] == (exit_code == 124)
    assert json.loads((directory / "frontier_report.json").read_text()) == result
    assert (directory / "release.json").exists()


def test_read_point_encloses_frontier_and_rechecks_membership(monkeypatch, evidence):
    points, done = evidence
    initial = deepcopy(points[0])
    extend(points, done)
    expected = points[0]
    monkeypatch.setattr(module, "read_hierarchy", lambda *args: deepcopy(initial))
    monkeypatch.setattr(module, "frontier_snapshot", lambda root, target: deepcopy(expected["frontier"]))
    monkeypatch.setattr(module, "host_snapshot", lambda: deepcopy(expected["host"][1]))
    original_read = module.Path.read_text
    monkeypatch.setattr(module.Path, "read_text", lambda path, *a, **kw:
                        initial["native_membership"] if str(path) == "/proc/123/cgroup"
                        else original_read(path, *a, **kw))
    assert module.read_frontier_point(123, initial["native_membership"], 21816) == expected
    with pytest.raises(ValueError, match="disappeared or changed"):
        module.read_frontier_point(123, "different membership", 21816)


def test_frontier_failure_releases_worker_without_admitting_result(tmp_path, monkeypatch, evidence):
    import tests.unit.test_measure_native_hierarchy_step as tests
    monkeypatch.setattr(tests, "module", module)
    directory = setup_measure(tmp_path, monkeypatch, evidence)
    def fail(*args):
        assert args[-1] == directory / "failed_point.json"
        raise ValueError("Frontier topology changed")
    monkeypatch.setattr(module, "read_frontier_point", fail)
    with pytest.raises(ValueError, match="topology changed"):
        module.measure(["/usr/bin/true"], directory, 21816, 20, 96*1024**3, 60, 1.)
    assert (directory / "release.json").exists()
    assert (directory / "go.json").exists()
    assert not (directory / "frontier_report.json").exists()


def test_snapshot_failure_saved_before_reraise(tmp_path, monkeypatch, evidence):
    points, done = evidence
    monkeypatch.setattr(module, "read_hierarchy", lambda *args: deepcopy(points[0]))
    failure = module.FrontierSnapshotError(ValueError("changed"), {"inventory_before": {}, "inventory_after": {}})
    def fail(*args):
        raise failure
    monkeypatch.setattr(module, "frontier_snapshot", fail)
    path = tmp_path / "failed_point.json"
    with pytest.raises(module.FrontierSnapshotError, match="changed"):
        module.read_frontier_point(123, "membership", 21816, path)
    saved = json.loads(path.read_text())
    assert saved["frontier"] == failure.evidence
    assert saved["scientific_timings_admitted"] is False
    assert saved["status"] == "invalid_frontier_observation"
    with pytest.raises(FileExistsError):
        module.read_frontier_point(123, "membership", 21816, path)

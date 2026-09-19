import json

import pytest

from benchmark_tools import measure_native_dual_bracket_step as module
from benchmark_tools.measure_native_frontier_step import evaluate as original_evaluate
from tests.unit.test_measure_native_hierarchy_step import evidence, setup_measure
from tests.unit.test_measure_native_frontier_step import extend
from tests.unit.test_frontier_pressure_integration import add_pressure


def prepare(evidence):
    points, done = evidence
    extend(points, done)
    add_pressure(points)
    return points, done


def test_original_screening_preserved_exactly(evidence):
    points, done = prepare(evidence)
    original = original_evaluate(points, done, 21816)
    result = module.evaluate(points, done, 21816)
    assert result["original_screening"] == original
    assert len(result["narrow_intervals"]) == len(points) - 1
    assert result["narrow_flagged_intervals"] == [i for i,r in enumerate(result["narrow_intervals"]) if not r["screen_passed"]]
    assert not result["scientific_timings_admitted"]


@pytest.mark.parametrize("code", [0, 7, 124])
def test_complete_worker_lifecycle_and_native_status(tmp_path, monkeypatch, evidence, code):
    import tests.unit.test_measure_native_hierarchy_step as helpers
    points, done = prepare(evidence)
    monkeypatch.setattr(helpers, "module", module)
    monkeypatch.setattr(module, "read_hierarchy", None, raising=False)
    directory = setup_measure(tmp_path, monkeypatch, evidence, code)
    pending = iter(points)
    monkeypatch.setattr(module, "read_point", lambda *a: next(pending))
    result = module.measure(["/usr/bin/true"], directory, 21816, 20, 96*1024**3, 60, 1.)
    assert result["native"]["exit_code"] == code
    assert result["native"]["timed_out"] == (code == 124)
    assert result["status"] == ("command_exited_zero" if code == 0 else "command_failed")
    assert json.loads((directory / "dual_bracket_report.json").read_text()) == result
    assert (directory / "release.json").exists()
    assert not result["scientific_timings_admitted"]


def test_observation_failure_releases_worker(tmp_path, monkeypatch, evidence):
    import tests.unit.test_measure_native_hierarchy_step as helpers
    prepare(evidence)
    monkeypatch.setattr(helpers, "module", module)
    monkeypatch.setattr(module, "read_hierarchy", None, raising=False)
    directory = setup_measure(tmp_path, monkeypatch, evidence)
    def fail(*a):
        raise ValueError("identity changed")
    monkeypatch.setattr(module, "read_point", fail)
    with pytest.raises(ValueError, match="identity changed"):
        module.measure(["/usr/bin/true"], directory, 21816, 20, 96*1024**3, 60, 1.)
    assert (directory / "release.json").exists()
    assert (directory / "go.json").exists()
    assert not (directory / "dual_bracket_report.json").exists()


@pytest.mark.parametrize("fault", ["enclosure", "pressure", "frontier"])
def test_original_failures_still_rejected(evidence, fault):
    points, done = prepare(evidence)
    if fault == "enclosure":
        done["started_ns"] = 0
    elif fault == "pressure":
        del points[1]["native_pressure"]
    else:
        for key in ("inventory_before", "inventory_after"):
            points[1]["frontier"][key]["identities"]["/user.slice"][1] += 1
    with pytest.raises(ValueError):
        module.evaluate(points, done, 21816)

from copy import deepcopy
import inspect
import json

import pytest

from benchmark_tools import measure_lineage_boundary_step as module
from benchmark_tools import measure_native_lineage_step as periodic
from tests.unit.test_measure_native_hierarchy_step import evidence, setup_measure
from tests.unit.test_measure_native_lineage_step import extend_lineage


def test_same_worker_reader_and_polling_except_periodic_reads():
    expected = inspect.getsource(periodic.measure).replace(
        'points.append(read_point(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json"))\n'
        '                save(directory / f"point_{index:04d}.json", points[-1])\n'
        '                if completed:\n                    break',
        'if completed:\n'
        '                    points.append(read_point(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json"))\n'
        '                    save(directory / f"point_{index:04d}.json", points[-1])\n'
        '                    break')
    expected = expected.replace("lineage_report.json", "lineage_boundary_report.json")
    expected = expected.replace("Complete-command lineage engineering measurement, not controlled comparative timing.",
        "Boundary-only lineage engineering measurement, not controlled comparative timing.")
    expected = expected.replace("Original screens and failures retained; no threshold changes or wall-time correction.",
        "No periodic interval screening; missing interval flags do not establish quietness.")
    assert inspect.getsource(module.measure) == expected
    assert module.worker is periodic.worker and module.read_point is periodic.read_point


def test_same_whole_screen_without_interval_claim(evidence):
    points, done = extend_lineage(evidence)
    result = module.evaluate(points, done, 21816)
    expected = periodic.evaluate(points, done, 21816)
    assert result["whole_command_screen"] == expected["original_screening"]["original_threshold_screen"]["whole_command_screen"]
    assert result["lineage_boundary"] == expected["observation_window"]["lineage"]
    assert result["native_pressure_whole_command"] == expected["observation_window"]["native_pressure"]
    assert result["flagged_intervals"] is None and not result["interval_screening_available"]
    assert not result["scientific_timings_admitted"]


@pytest.mark.parametrize("fault", ["missing", "extra", "start", "finish", "identity", "pressure", "schema"])
def test_invalid_boundary_evidence(evidence, fault):
    points, done = extend_lineage(evidence)
    if fault == "missing":
        points.pop()
    elif fault == "extra":
        points.append(deepcopy(points[-1]))
    elif fault == "start":
        done["started_ns"] = points[0]["host"][1]["finished_monotonic_ns"]
    elif fault == "finish":
        done["finished_ns"] = points[-1]["host"][0]["started_monotonic_ns"]
    elif fault == "identity":
        for key in ("identities_before", "identities_after"):
            points[-1]["lineage"][key]["/"][1] += 1
    elif fault == "pressure":
        del points[-1]["native_pressure"]
    else:
        points[-1]["schema"] = "wrong"
    with pytest.raises(ValueError):
        module.evaluate(points, done, 21816)


@pytest.mark.parametrize("exit_code", [0, 7, 124])
def test_only_two_reads_despite_multiple_polls(tmp_path, monkeypatch, evidence, exit_code):
    import tests.unit.test_measure_native_hierarchy_step as tests
    points, done = extend_lineage(evidence)
    monkeypatch.setattr(tests, "module", module)
    monkeypatch.setattr(module, "read_hierarchy", None, raising=False)
    directory = setup_measure(tmp_path, monkeypatch, evidence, exit_code)
    reads, polls = [], []
    def read(*args):
        reads.append(args)
        return points[len(reads) - 1]
    def sleep(seconds):
        polls.append(seconds)
        if len(polls) == 3:
            (directory / "done.json").write_text(json.dumps(done))
    monkeypatch.setattr(module, "read_point", read)
    monkeypatch.setattr(module.time, "sleep", sleep)
    result = module.measure(["/usr/bin/true"], directory, 21816, 20, 96 * 1024**3, 60, 1.)
    assert len(reads) == 2 and len(polls) == 3
    assert sorted(p.name for p in directory.glob("point_*.json")) == ["point_0000.json", "point_0003.json"]
    assert result["native"]["exit_code"] == exit_code
    assert result["native"]["timed_out"] == (exit_code == 124)
    assert (directory / "release.json").exists()
    assert json.loads((directory / "lineage_boundary_report.json").read_text()) == result


def test_failed_read_releases_worker(tmp_path, monkeypatch, evidence):
    import tests.unit.test_measure_native_hierarchy_step as tests
    monkeypatch.setattr(tests, "module", module)
    monkeypatch.setattr(module, "read_hierarchy", None, raising=False)
    directory = setup_measure(tmp_path, monkeypatch, evidence)
    def fail(*args):
        raise ValueError("read failed")
    monkeypatch.setattr(module, "read_point", fail)
    with pytest.raises(ValueError):
        module.measure(["/usr/bin/true"], directory, 21816, 20, 96 * 1024**3, 60, 1.)
    assert (directory / "release.json").exists()
    assert not (directory / "lineage_boundary_report.json").exists()

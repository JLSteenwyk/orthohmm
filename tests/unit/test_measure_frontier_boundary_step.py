from copy import deepcopy
import inspect
import json
from pathlib import Path

import pytest

from benchmark_tools import measure_frontier_boundary_step as module
from benchmark_tools import measure_native_frontier_step as original
from tests.unit.test_measure_native_hierarchy_step import evidence, setup_measure
from tests.unit.test_measure_native_frontier_step import extend


def test_worker_and_cleanup_unchanged_except_periodic_reads():
    expected = inspect.getsource(original.measure)
    expected = expected.replace(
        'points.append(reader(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json"))\n'
        '                save(directory / f"point_{index:04d}.json", points[-1])\n'
        '                if completed:\n                    break',
        'if completed:\n'
        '                    points.append(reader(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json"))\n'
        '                    save(directory / f"point_{index:04d}.json", points[-1])\n'
        '                    break')
    expected = expected.replace("frontier_report.json", "boundary_report.json")
    expected = expected.replace("Complete-command frontier engineering test", "Boundary-only incremental-overhead control")
    expected = expected.replace("Original CPU thresholds retained; hierarchy residuals do not alter flags or correct wall time.",
                                "No periodic interval screening is available; absence of flags is not a quiet-host result.")
    assert inspect.getsource(module.measure) == expected
    assert module.worker is original.worker


def test_same_whole_command_screen_without_claiming_interval_coverage(evidence):
    points, done = evidence
    extend(points, done)
    expected = original.evaluate(points, done, 21816)
    result = module.evaluate(points, done, 21816)
    assert result["whole_command_screen"] == expected["original_threshold_screen"]["whole_command_screen"]
    assert result["hierarchy_boundary"] == expected["hierarchy_intervals"][0]
    assert result["frontier_boundary"] == expected["frontier_intervals"][0]
    assert result["flagged_intervals"] is None
    assert result["interval_screening_available"] is False
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("index", [0, 1, 2])
def test_longer_retained_boundaries_replay_but_are_not_off_arm_measurements(index):
    root = Path(__file__).resolve().parents[2]
    report = json.loads((root / "benchmark_tools/results/dgx_frontier_native_smokes_21831.json").read_text())
    measured = report["runs"][index]["verification"]["measurement"]
    points = [measured["points"][0], measured["points"][-1]]
    result = module.evaluate(points, measured["native"], measured["job_id"])
    assert result["whole_command_screen"] == measured["screening"]["original_threshold_screen"]["whole_command_screen"]
    assert result["interval_screening_available"] is False
    assert result["flagged_intervals"] is None
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["missing", "extra", "start", "finish", "boot", "target"])
def test_invalid_boundary_evidence_rejected(evidence, fault):
    points, done = evidence
    extend(points, done)
    if fault == "missing":
        points.pop()
    elif fault == "extra":
        points.append(deepcopy(points[-1]))
    elif fault == "start":
        done["started_ns"] = points[0]["host"][1]["finished_monotonic_ns"]
    elif fault == "finish":
        done["finished_ns"] = points[-1]["host"][0]["started_monotonic_ns"]
    elif fault == "boot":
        points[-1]["frontier"]["boot_id"] = "wrong"
    else:
        points[-1]["frontier"]["target"] = "/user.slice"
    with pytest.raises(ValueError):
        module.evaluate(points, done, 21816)


@pytest.mark.parametrize("exit_code", [0, 7, 124])
def test_boundary_only_reads_despite_multiple_completion_polls(tmp_path, monkeypatch, evidence, exit_code):
    import tests.unit.test_measure_native_hierarchy_step as tests

    points, done = evidence
    extend(points, done)
    monkeypatch.setattr(tests, "module", module)
    monkeypatch.setattr(module, "read_hierarchy", None, raising=False)
    directory = setup_measure(tmp_path, monkeypatch, evidence, exit_code)
    reads, polls = [], []
    def read(*args):
        reads.append(args)
        return points[len(reads)-1]
    def sleep(seconds):
        polls.append(seconds)
        if len(polls) == 3:
            (directory / "done.json").write_text(json.dumps(done))
    monkeypatch.setattr(module, "read_frontier_point", read)
    monkeypatch.setattr(module.time, "sleep", sleep)
    result = module.measure(["/usr/bin/true"], directory, 21816, 20, 96*1024**3, 60, 1.)
    assert len(reads) == 2
    assert len(polls) == 3
    assert sorted(p.name for p in directory.glob("point_*.json")) == ["point_0000.json", "point_0003.json"]
    assert result["native"]["exit_code"] == exit_code
    assert result["native"]["timed_out"] == (exit_code == 124)
    assert result["status"] == ("command_exited_zero" if exit_code == 0 else "command_failed")
    assert json.loads((directory / "boundary_report.json").read_text()) == result
    assert (directory / "release.json").exists()


def test_boundary_read_failure_releases_owned_worker(tmp_path, monkeypatch, evidence):
    import tests.unit.test_measure_native_hierarchy_step as tests

    monkeypatch.setattr(tests, "module", module)
    monkeypatch.setattr(module, "read_hierarchy", None, raising=False)
    directory = setup_measure(tmp_path, monkeypatch, evidence)
    def fail(*args):
        assert args[-1] == directory / "failed_point.json"
        raise ValueError("frontier changed")
    monkeypatch.setattr(module, "read_frontier_point", fail)
    with pytest.raises(ValueError, match="frontier changed"):
        module.measure(["/usr/bin/true"], directory, 21816, 20, 96*1024**3, 60, 1.)
    assert (directory / "go.json").exists()
    assert (directory / "release.json").exists()
    assert not (directory / "boundary_report.json").exists()

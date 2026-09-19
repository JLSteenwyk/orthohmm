import json

import pytest

from benchmark_tools import replay_dual_native_measurement as module
from tests.unit.test_measure_native_dual_bracket_step import prepare
from tests.unit.test_measure_native_hierarchy_step import evidence


def save(path, value):
    path.write_text(json.dumps(value))


@pytest.fixture
def archive(tmp_path, evidence):
    points, done = prepare(evidence)
    memory = dict(scope=module.interval_point(points[-1], 21816)["native_cpu_scope"],
                  errors=[], raw={"memory.current": "100", "memory.peak": "200"},
                  started_ns=done["finished_ns"] + 10, finished_ns=done["finished_ns"] + 20)
    measured = dict(job_id=21816, native=done, step_memory=memory, points=points,
                    scientific_timings_admitted=False, controlled_workload_verified=False,
                    publication_ready=False, status="command_exited_zero",
                    native_wall_s=(done["finished_ns"] - done["started_ns"]) / 1e9,
                    screening=module.evaluate(points, done, 21816))
    save(tmp_path / "command.json", dict(command=["/usr/bin/true"], cpus=20,
                                         timeout_s=900, interval_s=1.))
    return tmp_path, measured


def write(archive):
    directory, measured = archive
    save(directory / "dual_bracket_report.json", measured)
    save(directory / "done.json", measured["native"])
    save(directory / "step_memory.json", measured["step_memory"])
    for index, point in enumerate(measured["points"]):
        save(directory / f"point_{index:04d}.json", point)


def replay(archive):
    return module.replay(archive[0], 21816, ["/usr/bin/true"])


def test_replays_both_screens_without_admitting_timing(archive):
    write(archive)
    result = replay(archive)
    assert result["screening"] == archive[1]["screening"]
    assert result["original_flagged_intervals"] == result["screening"]["original_screening"]["original_threshold_screen"]["flagged_intervals"]
    assert result["narrow_flagged_intervals"] == result["screening"]["narrow_flagged_intervals"]
    assert len(result["evidence"]) == 6
    assert not result["scientific_timings_admitted"]


@pytest.mark.parametrize("fault", ["command", "cadence", "job", "exit", "timeout", "wall",
                                  "memory_scope", "memory_time", "memory_error", "memory_peak",
                                  "memory_fraction", "admission", "outer", "narrow", "pressure",
                                  "frontier", "point_gap", "point_extra", "raw_disagreement"])
def test_rejects_inconsistent_or_invalid_evidence(archive, fault):
    directory, measured = archive
    if fault == "job":
        measured["job_id"] += 1
    elif fault == "exit":
        measured["native"]["exit_code"] = 7
    elif fault == "timeout":
        measured["native"]["timed_out"] = True
    elif fault == "wall":
        measured["native_wall_s"] += 1
    elif fault.startswith("memory_"):
        memory = measured["step_memory"]
        if fault == "memory_scope":
            memory["scope"] += "/wrong"
        elif fault == "memory_time":
            memory["started_ns"] = measured["native"]["finished_ns"]
        elif fault == "memory_error":
            memory["errors"] = ["read failed"]
        else:
            memory["raw"]["memory.peak"] = "99" if fault == "memory_peak" else "200.1"
    elif fault == "admission":
        measured["publication_ready"] = True
    elif fault == "outer":
        measured["screening"]["original_screening"]["original_threshold_screen"]["flagged_intervals"].append(99)
    elif fault == "narrow":
        measured["screening"]["narrow_flagged_intervals"].append(99)
    elif fault == "pressure":
        del measured["points"][0]["native_pressure"]
    elif fault == "frontier":
        for key in ("inventory_before", "inventory_after"):
            measured["points"][1]["frontier"][key]["identities"]["/user.slice"][1] += 1
    write(archive)
    if fault in ("command", "cadence"):
        path = directory / "command.json"
        value = json.loads(path.read_text())
        value["command" if fault == "command" else "interval_s"] = ["/usr/bin/false"] if fault == "command" else True
        save(path, value)
    elif fault == "point_gap":
        (directory / "point_0001.json").rename(directory / "point_0002.json")
    elif fault == "point_extra":
        save(directory / "point_extra.json", {})
    elif fault == "raw_disagreement":
        save(directory / "done.json", {})
    with pytest.raises(ValueError):
        replay(archive)


def test_detects_evidence_changed_during_replay(archive, monkeypatch):
    write(archive)
    evaluate = module.evaluate
    def change(*args):
        result = evaluate(*args)
        save(archive[0] / "done.json", {})
        return result
    monkeypatch.setattr(module, "evaluate", change)
    with pytest.raises(ValueError):
        replay(archive)


def test_detects_point_added_during_replay(archive, monkeypatch):
    write(archive)
    evaluate = module.evaluate
    def change(*args):
        result = evaluate(*args)
        save(archive[0] / "point_0002.json", {})
        return result
    monkeypatch.setattr(module, "evaluate", change)
    with pytest.raises(ValueError, match="inventory changed"):
        replay(archive)


def test_explicit_control_timeout_preserves_original_default(archive):
    write(archive)
    path = archive[0] / "command.json"
    value = json.loads(path.read_text())
    value["timeout_s"] = 60
    save(path, value)
    with pytest.raises(ValueError, match="command, resources"):
        replay(archive)
    result = module.replay(archive[0], 21816, ["/usr/bin/true"], expected_timeout_s=60)
    assert result["status"] == "dual_native_measurement_replayed"


@pytest.mark.parametrize("timeout", [True, 60., 61, 0])
def test_unknown_timeout_rejected(tmp_path, timeout):
    with pytest.raises(ValueError, match="frozen diagnostic timeout"):
        module.replay(tmp_path, 1, [], expected_timeout_s=timeout)

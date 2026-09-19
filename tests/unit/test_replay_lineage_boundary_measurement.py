import json

import pytest

from benchmark_tools import replay_lineage_boundary_measurement as module
from tests.unit.test_measure_native_hierarchy_step import evidence
from tests.unit.test_measure_native_lineage_step import extend_lineage


def save(path, value):
    path.write_text(json.dumps(value))


@pytest.fixture
def archive(tmp_path, evidence):
    points, done = extend_lineage(evidence)
    memory = dict(scope=module.interval_point(points[-1], 21816)["native_cpu_scope"],
        errors=[], raw={"memory.current": "100", "memory.peak": "200"},
        started_ns=done["finished_ns"] + 10, finished_ns=done["finished_ns"] + 20)
    measured = dict(job_id=21816, native=done, step_memory=memory, points=points,
        scientific_timings_admitted=False, controlled_workload_verified=False,
        publication_ready=False, status="command_exited_zero",
        native_wall_s=(done["finished_ns"] - done["started_ns"]) / 1e9,
        screening=module.evaluate(points, done, 21816))
    for name, value in (("lineage_boundary_report.json", measured), ("done.json", done),
                        ("step_memory.json", memory), ("point_0000.json", points[0]),
                        ("point_0003.json", points[1])):
        save(tmp_path / name, value)
    save(tmp_path / "command.json", dict(command=["/usr/bin/true"], cpus=20, timeout_s=900, interval_s=1.))
    return tmp_path


def replay(path):
    return module.replay(path, 21816, ["/usr/bin/true"])


def test_noncontiguous_poll_indices_are_valid_boundaries(archive):
    result = replay(archive)
    assert result["status"] == "lineage_boundary_measurement_replayed"
    assert result["flagged_intervals"] is None
    assert not result["interval_screening_available"]
    assert not result["scientific_timings_admitted"]


@pytest.mark.parametrize("fault", ["missing", "extra", "no_zero", "bad_suffix", "noncanonical", "command",
    "memory_peak", "raw_mismatch", "screen", "admission", "pressure", "schema"])
def test_invalid_boundary_archive(archive, fault):
    if fault == "missing":
        (archive / "point_0003.json").unlink()
    elif fault == "extra":
        save(archive / "point_0002.json", {})
    elif fault == "no_zero":
        (archive / "point_0000.json").rename(archive / "point_0001.json")
    elif fault in ("bad_suffix", "noncanonical"):
        (archive / "point_0003.json").rename(archive / ("point_x.json" if fault == "bad_suffix" else "point_3.json"))
    elif fault == "command":
        save(archive / "command.json", {})
    elif fault == "raw_mismatch":
        save(archive / "done.json", {})
    else:
        path = archive / "lineage_boundary_report.json"
        value = json.loads(path.read_text())
        if fault == "memory_peak":
            value["step_memory"]["raw"]["memory.peak"] = "99"
            save(archive / "step_memory.json", value["step_memory"])
        elif fault == "screen":
            value["screening"]["flagged_intervals"] = []
        elif fault == "admission":
            value["publication_ready"] = True
        else:
            if fault == "pressure":
                del value["points"][1]["native_pressure"]
            else:
                value["points"][1]["schema"] = "frontier"
            save(archive / "point_0003.json", value["points"][1])
        save(path, value)
    with pytest.raises((ValueError, KeyError)):
        replay(archive)


def test_modified_evidence_during_replay(archive, monkeypatch):
    original = module.evaluate
    def change(*args):
        result = original(*args)
        save(archive / "done.json", {})
        return result
    monkeypatch.setattr(module, "evaluate", change)
    with pytest.raises(ValueError):
        replay(archive)

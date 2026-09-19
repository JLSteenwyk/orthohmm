import json

import pytest

from benchmark_tools import replay_frontier_overhead_measurement as module
from tests.unit.test_frontier_pressure_integration import add_pressure
from tests.unit.test_measure_native_frontier_step import extend
from tests.unit.test_measure_native_hierarchy_step import evidence


def save(path, value):
    path.write_text(json.dumps(value))


@pytest.fixture(params=["boundary", "periodic"])
def archive(tmp_path, evidence, request):
    points, done = evidence
    extend(points, done)
    add_pressure(points)
    memory = dict(scope=module.interval_point(points[-1], 21816)["native_cpu_scope"],
                  errors=[], raw={"memory.current": "100", "memory.peak": "200"},
                  started_ns=done["finished_ns"] + 10, finished_ns=done["finished_ns"] + 20)
    mode = request.param
    evaluate = module.boundary_evaluate if mode == "boundary" else module.periodic_evaluate
    measured = dict(job_id=21816, native=done, step_memory=memory, points=points,
                    scientific_timings_admitted=False, controlled_workload_verified=False,
                    publication_ready=False, status="command_exited_zero",
                    native_wall_s=(done["finished_ns"] - done["started_ns"]) / 1e9,
                    screening=evaluate(points, done, 21816))
    save(tmp_path / "command.json", dict(command=["/usr/bin/true"], cpus=20,
                                         timeout_s=900, interval_s=1.))
    save(tmp_path / "done.json", done)
    save(tmp_path / "step_memory.json", memory)
    return tmp_path, mode, measured


def write_report(archive):
    directory, mode, measured = archive
    for index, point in enumerate(measured["points"]):
        save(directory / f"point_{index:04d}.json", point)
    name = "boundary_report.json" if mode == "boundary" else "frontier_report.json"
    save(directory / name, measured)


def replay(archive, expected=True):
    directory, mode, _ = archive
    write_report(archive)
    return module.replay(directory, mode, 21816, ["/usr/bin/true"],
                         expected_native_pressure=expected)


def test_pressure_enabled_replay_recomputes_diagnostics(archive):
    result = replay(archive)
    pressure = result["screening"]["native_pressure_whole_command"]
    assert pressure["native_stall_usec"]["cpu"] == {"some": 100, "full": 40}
    assert not result["scientific_timings_admitted"]
    if archive[1] == "boundary":
        assert result["flagged_intervals"] is None
        assert "native_pressure_intervals" not in result["screening"]
    else:
        assert len(result["screening"]["native_pressure_intervals"]) == 1


@pytest.mark.parametrize("all_points", [False, True])
def test_missing_pressure_cannot_be_downgraded_to_legacy(archive, all_points):
    _, mode, measured = archive
    points = measured["points"] if all_points else measured["points"][:1]
    for point in points:
        del point["native_pressure"]
    if all_points:
        evaluate = module.boundary_evaluate if mode == "boundary" else module.periodic_evaluate
        measured["screening"] = evaluate(measured["points"], measured["native"], 21816)
    with pytest.raises(ValueError, match="pressure coverage"):
        replay(archive)
    if all_points:
        assert replay(archive, False)["scientific_timings_admitted"] is False


def test_unexpected_pressure_rejected(archive):
    with pytest.raises(ValueError, match="pressure coverage"):
        replay(archive, False)


@pytest.mark.parametrize("expected", [0, 1, "true", []])
def test_nonboolean_expectation_rejected(archive, expected):
    with pytest.raises(ValueError, match="boolean"):
        replay(archive, expected)


@pytest.mark.parametrize("fault", ["reported_delta", "raw_totals", "identity", "scope"])
def test_pressure_replay_rejects_inconsistent_evidence(archive, fault):
    measured = archive[2]
    pressure = measured["points"][-1]["native_pressure"]
    if fault == "reported_delta":
        measured["screening"]["native_pressure_whole_command"]["native_stall_usec"]["cpu"]["some"] += 1
    elif fault == "raw_totals":
        pressure["pressure"]["cpu"]["totals"]["some"] += 1
    elif fault == "identity":
        pressure["scope_identity"][1] += 1
    else:
        pressure["scope"] += "/wrong"
    with pytest.raises(ValueError):
        replay(archive)

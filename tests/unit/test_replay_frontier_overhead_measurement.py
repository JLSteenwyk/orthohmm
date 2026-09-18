import copy
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import replay_frontier_overhead_measurement as module


@pytest.fixture(params=["raw_archive", "portable_report_fixture"])
def native(tmp_path, request):
    root = Path(__file__).resolve().parents[2]
    source = root / "benchmarks/work/dgx_frontier_native_21831/frontier_native_smoke_v1/run_00/measurement"
    directory = tmp_path / "measurement"
    if request.param == "raw_archive":
        if not source.exists():
            pytest.skip("Completed engineering raw archive is not present")
        shutil.copytree(source, directory)
        measured = json.loads((directory / "frontier_report.json").read_text())
        command = json.loads((directory / "command.json").read_text())["command"]
    else:
        from benchmark_tools.gnu_time_companion import command as timed_command
        from benchmark_tools.run_dgx_frontier_native_smoke import relocate
        results = root / "benchmark_tools/results"
        report = json.loads((results / "dgx_frontier_native_smokes_21831.json").read_text())
        measured = report["runs"][0]["verification"]["measurement"]
        spec = json.loads((results / "dgx_launcher_smoke_spec_20260917.json").read_text())
        run = relocate(spec["runs"][0])
        command = timed_command(run["native_argv"], str(Path(run["measurement_directory"]).parent / "native.time.tsv"))
        directory.mkdir()
        save(directory / "frontier_report.json", measured)
        save(directory / "command.json", dict(command=command, cpus=20, timeout_s=900, interval_s=1.))
        save(directory / "done.json", measured["native"])
        save(directory / "step_memory.json", measured["step_memory"])
        for index, point in enumerate(measured["points"]):
            save(directory / f"point_{index:04d}.json", point)
    return directory, measured, command


def save(path, value):
    path.write_text(json.dumps(value))


def test_retained_native_periodic_raw_evidence_replays(native):
    directory, measured, command = native
    result = module.replay(directory, "periodic", measured["job_id"], command)
    assert result["native_wall_s"] == measured["native_wall_s"]
    assert result["flagged_intervals"] == [1]
    assert result["whole_command_screen_passed"] is True
    assert result["scientific_timings_admitted"] is False


def test_synthetic_boundary_view_never_claims_interval_coverage(native):
    directory, measured, command = native
    paths = sorted(directory.glob("point_*.json"))
    for path in paths[1:-1]:
        path.unlink()
    measured["points"] = [measured["points"][0], measured["points"][-1]]
    measured["screening"] = module.boundary_evaluate(measured["points"], measured["native"], measured["job_id"])
    save(directory / "boundary_report.json", measured)
    result = module.replay(directory, "boundary", measured["job_id"], command)
    assert result["flagged_intervals"] is None
    assert result["screening"]["interval_screening_available"] is False
    # A synthetic view of prior periodic evidence is not an observer-off experiment.
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["command", "job", "native", "point", "missing_point", "extra_point",
                                    "wall", "flag", "admission", "memory_scope", "memory_counter",
                                    "memory_time", "memory_error", "failed_native", "boolean_interval"])
def test_corrupted_raw_evidence_rejected(native, fault):
    directory, measured, command = native
    job = measured["job_id"]
    if fault == "command":
        command = ["/different"]
    elif fault == "job":
        job += 1
    elif fault == "native":
        done = copy.deepcopy(measured["native"])
        done["exit_code"] = 7
        save(directory / "done.json", done)
    elif fault == "point":
        point = json.loads((directory / "point_0000.json").read_text())
        point["ticks"] += 1
        save(directory / "point_0000.json", point)
    elif fault == "missing_point":
        (directory / "point_0000.json").unlink()
    elif fault == "extra_point":
        shutil.copyfile(directory / "point_0000.json", directory / "point_9999.json")
    elif fault == "wall":
        measured["native_wall_s"] += .1
    elif fault == "flag":
        measured["screening"]["original_threshold_screen"]["flagged_intervals"] = []
    elif fault == "admission":
        measured["scientific_timings_admitted"] = True
    elif fault == "boolean_interval":
        payload = json.loads((directory / "command.json").read_text())
        payload["interval_s"] = True
        save(directory / "command.json", payload)
    elif fault.startswith("memory"):
        memory = measured["step_memory"]
        if fault == "memory_scope":
            memory["scope"] += "/wrong"
        elif fault == "memory_counter":
            memory["raw"]["memory.current"] = str(int(memory["raw"]["memory.peak"]) + 1)
        elif fault == "memory_time":
            memory["started_ns"] = measured["native"]["finished_ns"]
        else:
            memory["errors"] = ["read error"]
        save(directory / "step_memory.json", memory)
    else:
        measured["native"]["exit_code"] = 7
        measured["status"] = "command_failed"
        save(directory / "done.json", measured["native"])
    save(directory / "frontier_report.json", measured)
    with pytest.raises(ValueError):
        module.replay(directory, "periodic", job, command)

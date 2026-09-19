"""Replay dual-bracket native evidence without admitting comparative timings."""

import json
import math
from pathlib import Path

from benchmark_tools.measure_native_dual_bracket_step import evaluate, interval_point
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def replay(directory, job_id, expected_command, expected_timeout_s=900):
    directory = Path(directory)
    if type(expected_timeout_s) is not int or expected_timeout_s not in (60, 900):
        raise ValueError("Require a frozen diagnostic timeout")
    if type(job_id) is not int or job_id <= 0:
        raise ValueError("Expected job identity must be a positive integer")
    files = [directory / name for name in (
        "dual_bracket_report.json", "command.json", "done.json", "step_memory.json")]
    point_files = sorted(directory.glob("point_*.json"))
    if len(point_files) < 2 or point_files != [
            directory / f"point_{index:04d}.json" for index in range(len(point_files))]:
        raise ValueError("Raw point inventory must be contiguous from zero")
    evidence = [record(path) for path in [*files, *point_files]]
    measured, command, done, memory = [json.loads(path.read_text()) for path in files]
    if (command != dict(command=expected_command, cpus=20, timeout_s=expected_timeout_s, interval_s=1.)
            or type(command["cpus"]) is not int or type(command["timeout_s"]) is not int
            or type(command["interval_s"]) not in (int, float)):
        raise ValueError("Measured command, resources or cadence differ")
    if type(measured["job_id"]) is not int or measured["job_id"] != job_id:
        raise ValueError("Native job identity differs")
    if measured["native"] != done or measured["step_memory"] != memory:
        raise ValueError("Embedded and raw native/memory evidence disagree")
    points = [json.loads(path.read_text()) for path in point_files]
    if points != measured["points"]:
        raise ValueError("Raw points differ from report")
    if any(measured[key] is not False for key in (
            "scientific_timings_admitted", "controlled_workload_verified", "publication_ready")):
        raise ValueError("Unexpected timing admission")
    if (type(done["exit_code"]) is not int or done["exit_code"] != 0
            or done["timed_out"] is not False or measured["status"] != "command_exited_zero"):
        raise ValueError("Native process did not finish successfully")
    if any(type(done[key]) is not int or done[key] <= 0 for key in ("started_ns", "finished_ns")):
        raise ValueError("Invalid native time boundaries")
    wall = (done["finished_ns"] - done["started_ns"]) / 1e9
    if (not math.isfinite(wall) or wall <= 0 or type(measured["native_wall_s"]) not in (int, float)
            or measured["native_wall_s"] != wall):
        raise ValueError("Native wall duration does not reproduce")
    screening = evaluate(points, done, job_id)
    if screening != measured["screening"]:
        raise ValueError("Dual CPU replay differs from retained screening")
    if (memory["errors"] or memory["scope"] != interval_point(points[-1], job_id)["native_cpu_scope"]
            or any(type(memory[key]) is not int for key in ("started_ns", "finished_ns"))
            or not done["finished_ns"] < memory["started_ns"] <= memory["finished_ns"]):
        raise ValueError("Final memory observation has wrong scope, time or read errors")
    counters = [memory["raw"][key].strip() for key in ("memory.current", "memory.peak")]
    if any(not value.isascii() or not value.isdecimal() for value in counters):
        raise ValueError("Invalid final memory counters")
    current, peak = map(int, counters)
    if current > peak:
        raise ValueError("Invalid final memory counters")
    for item in evidence:
        check(item)
    if point_files != sorted(directory.glob("point_*.json")):
        raise ValueError("Raw point inventory changed during replay")
    original = screening["original_screening"]["original_threshold_screen"]
    return dict(status="dual_native_measurement_replayed", native_wall_s=wall,
        original_flagged_intervals=original["flagged_intervals"],
        narrow_flagged_intervals=screening["narrow_flagged_intervals"],
        screening=screening, memory=memory, measured=measured, evidence=evidence,
        source=record(__file__), scientific_timings_admitted=False,
        controlled_workload_verified=False, publication_ready=False,
        limitations=["Measurement replay only; scheduler, authorized command, runtime and native outputs require separate audit.",
                     "Observation windows include wrapper work and are not exact command boundaries.",
                     "Final cgroup memory includes wrappers/cache and is not maximum process RSS.",
                     "Neither residual screen identifies foreign workloads or establishes controlled timing."])

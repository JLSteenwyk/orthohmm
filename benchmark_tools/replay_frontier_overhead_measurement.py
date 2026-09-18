"""Replay archived overhead measurements without making timing admission claims."""

import json
import math
from pathlib import Path

from benchmark_tools.measure_native_frontier_step import evaluate as periodic_evaluate, interval_point
from benchmark_tools.measure_frontier_boundary_step import evaluate as boundary_evaluate
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def replay(directory, mode, job_id, expected_command):
    if mode not in {"boundary", "periodic"}:
        raise ValueError("Unknown overhead observation mode")
    report_name = "boundary_report.json" if mode == "boundary" else "frontier_report.json"
    files = [directory / name for name in (report_name, "command.json", "done.json", "step_memory.json")]
    point_files = sorted(directory.glob("point_*.json"))
    evidence = [record(path) for path in [*files, *point_files]]
    measured, command, done, memory = [json.loads(path.read_text()) for path in files]
    if (command != dict(command=expected_command, cpus=20, timeout_s=900, interval_s=1.)
            or type(command["interval_s"]) not in (int, float)):
        raise ValueError("Measured command, resources or cadence differ")
    if measured["job_id"] != job_id or type(measured["job_id"]) is not int:
        raise ValueError("Native job identity differs")
    if measured["native"] != done or measured["step_memory"] != memory:
        raise ValueError("Embedded and raw native/memory evidence disagree")
    points = [json.loads(path.read_text()) for path in point_files]
    if points != measured["points"]:
        raise ValueError("Raw point inventory differs from report")
    if any(measured[key] is not False for key in (
            "scientific_timings_admitted", "controlled_workload_verified", "publication_ready")):
        raise ValueError("Unexpected timing admission")
    if (type(done["exit_code"]) is not int or done["exit_code"] != 0 or done["timed_out"] is not False
            or measured["status"] != "command_exited_zero"):
        raise ValueError("Native process did not finish successfully")
    wall = (done["finished_ns"] - done["started_ns"]) / 1e9
    if (not math.isfinite(wall) or wall <= 0 or type(measured["native_wall_s"]) not in (int, float)
            or measured["native_wall_s"] != wall):
        raise ValueError("Native wall duration does not reproduce")
    evaluate = boundary_evaluate if mode == "boundary" else periodic_evaluate
    screening = evaluate(points, done, job_id)
    if screening != measured["screening"]:
        raise ValueError("CPU replay differs from retained screening")
    if (memory["errors"] or memory["scope"] != interval_point(points[-1], job_id)["native_cpu_scope"]
            or not done["finished_ns"] < memory["started_ns"] <= memory["finished_ns"]):
        raise ValueError("Final memory observation has wrong scope, time or read errors")
    current, peak = (int(memory["raw"][key]) for key in ("memory.current", "memory.peak"))
    if not 0 <= current <= peak:
        raise ValueError("Invalid final memory counters")
    for item in evidence:
        check(item)
    if point_files != sorted(directory.glob("point_*.json")):
        raise ValueError("Raw point inventory changed during replay")
    original = screening if mode == "boundary" else screening["original_threshold_screen"]
    return dict(status="overhead_measurement_replayed", native_wall_s=wall,
        whole_command_screen_passed=original["whole_command_screen"]["screen_passed"],
        flagged_intervals=original["flagged_intervals"], screening=screening,
        memory=memory, measured=measured, evidence=evidence, source=record(__file__),
        scientific_timings_admitted=False,
        limitations=["Measurement replay only; scheduler, command authorization, runtime and native outputs require separate audit.",
                     "Final cgroup memory includes wrappers/cache and is not maximum process RSS."])

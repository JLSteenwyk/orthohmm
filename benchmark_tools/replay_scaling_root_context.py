"""Replay long-run raw measurements, retaining native failures without admission."""

import json
import math
from pathlib import Path

from benchmark_tools.measure_native_lineage_step import evaluate as evaluate_lineage, interval_point
from benchmark_tools.measure_native_root_context import evaluate, lineage_identity
from benchmark_tools.measure_scaling_root_context import TIMEOUT, validate
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.slurm_resource_snapshot import counters
from benchmark_tools.verify_lineage_native_provenance import same
from benchmark_tools.replay_host_process_observation import replay as replay_host
from benchmark_tools.slurm_resource_snapshot import scoped_path


def native_outcome(measured, done):
    code, timed_out = done["exit_code"], done["timed_out"]
    if type(code) is not int or type(timed_out) is not bool or (timed_out and code != 124):
        raise ValueError("Invalid native exit/timeout combination")
    expected_status = "command_exited_zero" if code == 0 else "command_failed"
    if measured["status"] != expected_status:
        raise ValueError("Collector status contradicts native exit")
    if any(type(done[k]) is not int or done[k] <= 0 for k in ("started_ns", "finished_ns")):
        raise ValueError("Invalid native time boundaries")
    wall = (done["finished_ns"]-done["started_ns"])/1e9
    if (not 0 < wall <= TIMEOUT + 30 or (timed_out and wall < TIMEOUT)
            or type(measured["native_wall_s"]) not in (int, float)
            or not math.isfinite(measured["native_wall_s"]) or measured["native_wall_s"] != wall):
        raise ValueError("Native wall duration violates the recorded boundary")
    return ("timed_out" if timed_out else "exited_zero" if code == 0 else "exited_nonzero"), wall


def point_inventory(directory):
    paths = sorted(directory.glob("point_*.json"))
    if len(paths) < 2 or paths != [directory / f"point_{i:06d}.json" for i in range(len(paths))]:
        raise ValueError("Require contiguous six-digit point inventory from zero")
    return paths


def replay(directory, job_id, expected_command):
    directory = Path(directory).absolute()
    if type(job_id) is not int or job_id <= 0:
        raise ValueError("Require positive expected job identity")
    validate(expected_command, 20, TIMEOUT, 1.)
    entries = sorted(directory.iterdir())
    if any((directory / name).exists() or (directory / name).is_symlink() for name in (
            "aborted_before_native.json", "failed_point.json", "failed_root_context.json")):
        raise ValueError("Failure/abort marker contradicts complete measurement replay")
    files = [directory / name for name in (
        "lineage_report.json", "command.json", "done.json", "step_memory.json", "root_context_report.json",
        "go.json", "release.json")]
    point_files = point_inventory(directory)
    if directory.is_symlink() or any(p.is_symlink() for p in [*files, *point_files]):
        raise ValueError("Direct evidence symlinks are not supported")
    evidence = [record(p) for p in [*files, *point_files]]
    measured, command, done, memory, context, go, release = [json.loads(p.read_text()) for p in files]
    if not same(go, {"go": True}) or not same(release, {"release": True}):
        raise ValueError("Native handoff/release gates differ")
    if not same(command, dict(command=expected_command, cpus=20, timeout_s=TIMEOUT, interval_s=1.)):
        raise ValueError("Command, resources, timeout or cadence differ")
    if type(measured["job_id"]) is not int or measured["job_id"] != job_id:
        raise ValueError("Native job identity differs")
    if not same(measured["native"], done) or not same(measured["step_memory"], memory):
        raise ValueError("Embedded native/memory data disagree with raw evidence")
    if any(measured[k] is not False for k in ("scientific_timings_admitted", "controlled_workload_verified", "publication_ready")):
        raise ValueError("Unexpected timing admission")
    outcome, wall = native_outcome(measured, done)
    points = [json.loads(p.read_text()) for p in point_files]
    if not same(points, measured["points"]):
        raise ValueError("Raw points differ from report")
    host_paths = [directory / name for name in ("host_processes.jsonl", "host_process_summary.json")]
    host_replay = None
    if "host_process_observation" in measured or any(p.exists() or p.is_symlink() for p in host_paths):
        if any(p.is_symlink() or not p.is_file() for p in host_paths):
            raise ValueError("Incomplete or indirect host process evidence")
        evidence.extend(record(p) for p in host_paths)
        summary = json.loads(host_paths[1].read_text())
        if not same(summary, measured.get("host_process_observation")):
            raise ValueError("Embedded host process summary differs")
        scope = scoped_path(points[0]["native_membership"], job_id)
        job_scope = next(p for p in scope.parents if p.name == f"job_{job_id}")
        host_replay = replay_host(host_paths[0], summary, str(job_scope),
            done["started_ns"] / 1e9, done["finished_ns"] / 1e9)
    screening = evaluate_lineage(points, done, job_id)
    if not same(screening, measured["screening"]):
        raise ValueError("Lineage screening does not reproduce")
    if (memory["errors"] or memory["scope"] != interval_point(points[-1], job_id)["native_cpu_scope"]
            or any(type(memory[k]) is not int for k in ("started_ns", "finished_ns"))
            or not max(done["finished_ns"], points[-1]["root_context"]["host_after"]["finished_ns"])
                   < memory["started_ns"] <= memory["finished_ns"]):
        raise ValueError("Invalid final memory scope/read window")
    values = [memory["raw"][k].strip() for k in ("memory.current", "memory.peak")]
    if any(not v.isascii() or not v.isdecimal() for v in values) or int(values[0]) > int(values[1]):
        raise ValueError("Invalid memory gauge values")
    events = counters(memory["raw"]["memory.events"])
    if not {"low", "high", "max", "oom", "oom_kill"} <= events.keys():
        raise ValueError("Incomplete memory event counters")
    reproduced = dict(status="native_root_context_measured", job_id=job_id, native_wall_s=wall,
        context=evaluate(points, job_id), lineage_report=lineage_identity(directory),
        scientific_timings_admitted=False, environmental_validity_established=False)
    if not same(context, reproduced):
        raise ValueError("Supplementary context does not reproduce")
    for item in evidence:
        check(item)
    if entries != sorted(directory.iterdir()) or point_files != point_inventory(directory):
        raise ValueError("Raw point inventory changed during replay")
    return dict(status="scaling_root_context_measurement_replayed", native_outcome=outcome,
        native_exit_code=done["exit_code"], native_wall_s=wall, measured=measured, memory=memory,
        host_process_replay=host_replay,
        memory_events=events, screening=screening, context=reproduced["context"],
        original_flagged_intervals=screening["original_screening"]["original_threshold_screen"]["flagged_intervals"],
        narrow_flagged_intervals=screening["narrow_flagged_intervals"], evidence=evidence,
        source=record(__file__), scientific_timings_admitted=False, environmental_validity_established=False,
        native_outputs_validated=False, publication_ready=False,
        limitations=["Replay preserves successful, failed and timed-out native outcomes, not output validity.",
            "Separate authorization, scheduler/session, worker/recipe identity and environment audits remain mandatory.",
            "Raw observation windows and flags retained; no attribution, exclusion or overhead subtraction.",
            "Full raw points are loaded for replay; memory use and long-run execution still require assessment."])

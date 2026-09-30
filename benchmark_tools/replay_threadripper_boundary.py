"""Check boundary-control raw receipts; never validate outputs or admit timings."""

import json
from pathlib import Path

from benchmark_tools.disk_observation_sequence import DiskObservations
from benchmark_tools.measure_native_lineage_step import evaluate as evaluate_lineage, interval_point
from benchmark_tools.measure_native_root_context import evaluate as evaluate_context
from benchmark_tools.native_completion import completion_evidence
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_threadripper_allocation import validate as validate_allocation
from benchmark_tools.replay_host_process_observation import replay as replay_host
from benchmark_tools.replay_threadripper_scaling import (
    TIMEOUT, native_outcome, replay_affinity, point_inventory, validate_job_memory,
)
from benchmark_tools.report_finalization import validate as validate_finalization
from benchmark_tools.slurm_resource_snapshot import scoped_path, counters
from benchmark_tools.verify_lineage_native_provenance import same


def replay(directory, job_id, expected_command, *, expected_launcher, expected_worker,
           original_directory=None):
    directory = Path(directory).absolute()
    if directory.is_symlink() or directory.resolve() != directory or not directory.is_dir():
        raise ValueError("Require a direct measurement directory")
    if type(job_id) is not int or job_id <= 0:
        raise ValueError("Require positive expected job identity")
    from benchmark_tools.measure_threadripper_scaling import validate
    validate(expected_command, 32, TIMEOUT, 1.)
    for path in (expected_launcher, expected_worker):
        if not isinstance(path, str) or not Path(path).is_absolute() or ".." in Path(path).parts:
            raise ValueError("Require explicit absolute launcher and worker bindings")
    original = directory if original_directory is None else Path(original_directory)
    if not original.is_absolute() or ".." in original.parts:
        raise ValueError("Invalid original measurement directory")
    entries = sorted(directory.iterdir())
    markers = ("boundary_failure.json", "aborted_before_native.json", "failed_point.json",
               "failed_root_context.json", "release_freshness_failed.json")
    if any((directory / name).exists() or (directory / name).is_symlink() for name in markers):
        raise ValueError("Failure marker contradicts complete boundary replay")
    names = ("boundary_report.json", "command.json", "done.json", "step_memory.json",
             "go.json", "release.json", "ready.json", "native_completion.json",
             "job_memory_before.json", "job_memory_after.json", "host_processes.jsonl",
             "host_process_summary.json", "report_finalization.json")
    point_files = point_inventory(directory)
    if len(point_files) != 2:
        raise ValueError("Boundary control requires exactly two raw points")
    files = [directory / name for name in names] + point_files
    if any(p.is_symlink() or not p.is_file() for p in files):
        raise ValueError("Require complete direct raw boundary evidence")
    evidence = [record(p) for p in files]
    raw = {name: json.loads((directory / name).read_text()) for name in names if name.endswith(".json")}
    measured, command, done, memory = (raw[name] for name in (
        "boundary_report.json", "command.json", "done.json", "step_memory.json"))
    ready = raw["ready.json"]
    policy = dict(native_points=2, periodic_native_sampling=False,
                  completion_poll_interval_s=1., common_host_interval_s=30.)
    if (measured.get("schema") != "threadripper_boundary_control_v1"
            or measured.get("collector_arm") != "boundary" or not same(measured.get("policy"), policy)):
        raise ValueError("Unexpected boundary schema, arm or observation policy")
    if any(measured.get(k) is not False for k in ("scientific_timings_admitted",
            "controlled_workload_verified", "native_outputs_validated", "publication_ready")):
        raise ValueError("Unexpected boundary admission")
    if not same(raw["go.json"], {"go": True}) or not same(raw["release.json"], {"release": True}):
        raise ValueError("Native release gates differ")
    if not same(command, dict(command=expected_command, cpus=32, timeout_s=TIMEOUT, interval_s=1.)):
        raise ValueError("Native command or resource settings differ")
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1",
        "--cpus-per-task=64", "--cpu-bind=mask_cpu:0xffffffff", expected_launcher,
        "-B", expected_worker, "--worker", str(original)]
    if not same(measured.get("launched"), launched):
        raise ValueError("Worker launch differs from external bindings")
    if type(measured.get("job_id")) is not int or measured["job_id"] != job_id:
        raise ValueError("Boundary job identity differs")
    if not same(measured.get("native"), done) or not same(measured.get("step_memory"), memory):
        raise ValueError("Embedded native/memory data disagree")
    validate_allocation(ready["placement"], measured["placement"], 64)
    if (ready["placement"]["affinity"] != list(range(32))
            or ready["placement"]["slurm"]["SLURM_JOB_ID"] != str(job_id)):
        raise ValueError("Native placement or allocation job differs")
    outcome, wall = native_outcome(measured, done)
    points = DiskObservations(directory, 2)
    point_records = [dict(record(p), path=str(original / p.name)) for p in point_files]
    if "points" in measured or not same(measured.get("point_records"), point_records):
        raise ValueError("Raw point bindings differ")
    affinity = [replay_affinity(p, job_id) for p in points]
    if (points[0]["thread_affinity"]["finished_ns"] >= done["started_ns"]
            or points[-1]["host"][0]["started_monotonic_ns"] <= done["finished_ns"]):
        raise ValueError("Native boundary points do not enclose the command")
    if (points[0]["native_membership"] != ready["cgroup"]
            or points[1]["native_membership"] != ready["cgroup"]):
        raise ValueError("Native ready/point membership differs")
    completion = completion_evidence(points[0]["thread_affinity"], points[1]["thread_affinity"],
                                     ready["pid"], done["finished_ns"])
    if (completion["errors"] or not same(raw["native_completion.json"], completion)
            or not same(measured.get("native_completion"), completion)):
        raise ValueError("Native completion does not reproduce")
    screening, context = evaluate_lineage(points, done, job_id), evaluate_context(points, job_id)
    if not same(measured.get("screening"), screening) or not same(measured.get("root_context"), context):
        raise ValueError("Boundary screening/context does not reproduce")
    if (memory["errors"] or memory["scope"] != interval_point(points[-1], job_id)["native_cpu_scope"]
            or any(type(memory[k]) is not int for k in ("started_ns", "finished_ns"))
            or not max(done["finished_ns"], points[-1]["thread_affinity"]["finished_ns"])
                   < memory["started_ns"] <= memory["finished_ns"]):
        raise ValueError("Final native memory scope or window differs")
    gauges = [memory["raw"][k].strip() for k in ("memory.current", "memory.peak")]
    if (any(not v.isascii() or not v.isdecimal() for v in gauges)
            or int(gauges[0]) > int(gauges[1])):
        raise ValueError("Invalid native memory gauges")
    events = counters(memory["raw"]["memory.events"])
    if not {"low", "high", "max", "oom", "oom_kill"} <= events.keys():
        raise ValueError("Missing native memory event counters")
    before, after = raw["job_memory_before.json"], raw["job_memory_after.json"]
    if not same(measured.get("job_memory"), dict(before=before, after=after)):
        raise ValueError("Embedded job memory differs")
    scope = scoped_path(ready["cgroup"], job_id)
    job_scope = next(p for p in scope.parents if p.name == f"job_{job_id}")
    job_memory = validate_job_memory(before, after, job_scope, done, memory)
    finalization = raw["report_finalization.json"]
    final_memory = validate_finalization(finalization, job_id, after)
    validate_job_memory(before, final_memory, job_scope, done, memory)
    previous_events, final_events = (counters(r["raw"]["memory.events"]) for r in (after, final_memory))
    if (int(final_memory["raw"]["memory.peak"]) < int(after["raw"]["memory.peak"])
            or previous_events.keys() != final_events.keys()
            or any(final_events[k] < v for k, v in previous_events.items())):
        raise ValueError("Reporting-stage peak or memory event counters decreased")
    summary = raw["host_process_summary.json"]
    if not same(measured.get("host_process_observation"), summary):
        raise ValueError("Embedded host observation differs")
    host = replay_host(directory / "host_processes.jsonl", summary, str(job_scope),
                       done["started_ns"] / 1e9, done["finished_ns"] / 1e9)
    for item in evidence:
        check(item)
    if entries != sorted(directory.iterdir()) or point_files != point_inventory(directory):
        raise ValueError("Boundary inventory changed during replay")
    return dict(status="threadripper_boundary_measurement_replayed", native_outcome=outcome,
        native_exit_code=done["exit_code"], native_wall_s=wall, measured=measured,
        screening=screening, root_context=context, native_completion=completion,
        affinity_observation_statuses=affinity, memory_events=events, job_memory=job_memory,
        report_finalization=finalization, host_process_replay=host, evidence=evidence,
        source=record(__file__), scientific_timings_admitted=False,
        controlled_workload_verified=False, native_outputs_validated=False, publication_ready=False,
        limitations=["Replays two native boundary points, not continuous native containment or affinity.",
            "Successful replay is not zero native exit, output equality, runtime/source or scheduler admission.",
            "Whole-host observation is checked arithmetically, not a quiet-host or contamination-policy decision.",
            "Common host-monitor cost is retained; no overhead subtraction or total-cost estimate."])

"""Replay the new allocated-core schema without translating historical evidence."""

import json
from pathlib import Path

from benchmark_tools.measure_allocated_threadripper_scaling import SCHEMA, validate_command, selected_ids, sources
from benchmark_tools.native_factorial_allocated_placement import validate as validate_placement
from benchmark_tools.replay_threadripper_scaling import (
    replay_completion, native_outcome, point_inventory, validate_job_memory)
from benchmark_tools.measure_native_lineage_step import evaluate as evaluate_lineage, interval_point
from benchmark_tools.measure_native_root_context import evaluate, lineage_identity
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.slurm_resource_snapshot import counters, scoped_path
from benchmark_tools.verify_lineage_native_provenance import same
from benchmark_tools.replay_host_process_observation import replay as replay_host
from benchmark_tools.disk_observation_sequence import DiskObservations
from benchmark_tools.report_finalization import validate as validate_finalization


def replay_affinity(point, job, allowed):
    allowed = selected_ids(allowed)
    value = point["thread_affinity"]
    scope = scoped_path(point["native_membership"], job)
    while scope.name != "user":
        if scope.name.startswith("step_") or scope == scope.parent:
            raise ValueError("Missing user subtree")
        scope = scope.parent
    expected = str(Path("/sys/fs/cgroup") / str(scope).lstrip("/"))
    if (value["scope"] != expected or value["allowed_cpus"] != allowed
            or value["full_run_affinity_verified"] is not False
            or any(type(value[k]) is not int for k in ("started_ns", "finished_ns"))
            or not point["root_context"]["host_after"]["finished_ns"]
                   <= value["started_ns"] <= value["finished_ns"]):
        raise ValueError("Affinity scope, selected policy or window differs")
    for key in ("initial_tids", "final_tids"):
        ids = value[key]
        if (not isinstance(ids, list) or any(type(t) is not int or t <= 0 for t in ids)
                or ids != sorted(set(ids))):
            raise ValueError("Malformed thread inventory")
    seen, violations = [], []
    for row in value["threads"]:
        tid, cpus = row["tid"], row["affinity"]
        if (type(tid) is not int or tid not in value["initial_tids"] or tid in seen
                or type(row["start_ticks"]) is not int or row["start_ticks"] < 0
                or not cpus or any(type(c) is not int or c < 0 for c in cpus)
                or cpus != sorted(set(cpus))):
            raise ValueError("Malformed thread observation")
        if not scoped_path(row["cgroup"], job).is_relative_to(scope):
            raise ValueError("Observed thread outside native scope")
        outside = sorted(set(cpus) - set(allowed))
        if row["outside_cpus"] != outside:
            raise ValueError("Affinity violation does not reproduce")
        if outside:
            violations.append(tid)
        seen.append(tid)
    errors = value["errors"]
    if not isinstance(errors, list) or any(not isinstance(e, dict) or not e.get("error") for e in errors):
        raise ValueError("Malformed affinity errors")
    if (not value["initial_tids"] or seen != value["initial_tids"]
            or value["initial_tids"] != value["final_tids"]) and not errors:
        raise ValueError("Observation gap was hidden")
    status = "violation" if violations else "incomplete" if errors else "observed_within_affinity"
    if value["violating_tids"] != violations or value["status"] != status:
        raise ValueError("Affinity status does not reproduce")
    return status


def replay(directory, job_id, expected_command):
    directory = Path(directory).absolute()
    if type(job_id) is not int or job_id <= 0:
        raise ValueError("Require positive expected job identity")
    entries = sorted(directory.iterdir())
    if any((directory / name).exists() or (directory / name).is_symlink() for name in (
            "aborted_before_native.json", "failed_point.json", "failed_root_context.json",
            "release_freshness_failed.json")):
        raise ValueError("Failure/abort marker contradicts complete measurement replay")
    names = ("lineage_report.json", "command.json", "done.json", "step_memory.json",
        "root_context_report.json", "go.json", "release.json", "ready.json")
    files = [directory / name for name in names]
    point_files = point_inventory(directory)
    if directory.is_symlink() or any(p.is_symlink() or not p.is_file() for p in [*files, *point_files]):
        raise ValueError("Require direct measurement evidence")
    evidence = [record(p) for p in [*files, *point_files]]
    measured, command, done, memory, context, go, release, ready = [json.loads(p.read_text()) for p in files]
    if measured.get("schema") != SCHEMA:
        raise ValueError("Require allocation-aware collector schema; no historical translation")
    validate_command(command)
    if not same(command["command"], expected_command):
        raise ValueError("Expected native command differs")
    if not same(measured["sources"], sources()):
        raise ValueError("Collector sources differ")
    evidence.extend(measured["sources"])
    if not same(go, {"go": True}) or not same(release, {"release": True}):
        raise ValueError("Native handoff/release gates differ")
    if type(measured["job_id"]) is not int or measured["job_id"] != job_id:
        raise ValueError("Native job identity differs")
    if not same(measured["native"], done) or not same(measured["step_memory"], memory):
        raise ValueError("Embedded native/memory data disagree with raw evidence")
    if any(measured[k] is not False for k in (
            "scientific_timings_admitted", "controlled_workload_verified", "publication_ready")):
        raise ValueError("Unexpected timing admission")
    allowed = validate_placement(ready["allocated_placement"], job_id)
    if (not same(measured["allocated_placement"], ready["allocated_placement"])
            or not same(ready["placement"], ready["allocated_placement"]["bound"])
            or not same(measured["placement"], ready["placement"])
            or type(ready["pid"]) is not int or ready["pid"] != ready["placement"]["pid"]
            or ready["cgroup"] != ready["placement"]["cgroup"]):
        raise ValueError("Worker or placement binding differs")
    launched = measured["launched"]
    if (not isinstance(launched, list) or len(launched) != 12
            or launched[:7] != ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1",
                               "--cpus-per-task=64", "--cpu-bind=cores"]
            or not isinstance(launched[7], str) or not Path(launched[7]).is_absolute()
            or launched[8:] != ["-B", str(Path(__file__).with_name(
                "measure_allocated_threadripper_scaling.py").resolve()), "--worker", str(directory)]):
        raise ValueError("Allocated-core worker launch differs")
    outcome, wall = native_outcome(measured, done)
    points = DiskObservations(directory, len(point_files))
    if ("points" in measured or [points.path(i) for i in range(len(points))] != point_files
            or not same(points.records(), measured["point_records"])):
        raise ValueError("Disk-backed point identities differ")
    if ready["allocated_placement"]["finished_ns"] >= points[0]["host"][0]["started_monotonic_ns"]:
        raise ValueError("Placement does not precede observations")
    affinity = [replay_affinity(p, job_id, allowed) for p in points]
    completion, completion_record = replay_completion(directory, measured, points, ready, done)
    evidence.append(completion_record)
    for left, right in zip(points, points[1:]):
        if left["thread_affinity"]["finished_ns"] >= right["host"][0]["started_monotonic_ns"]:
            raise ValueError("Affinity observation overlaps the next point")
    scope = scoped_path(points[0]["native_membership"], job_id)
    if scope != scoped_path(ready["cgroup"], job_id):
        raise ValueError("Anchor membership differs from bound worker")
    job_scope = next(p for p in scope.parents if p.name == f"job_{job_id}")
    host_paths = [directory / name for name in ("host_processes.jsonl", "host_process_summary.json")]
    if any(p.is_symlink() or not p.is_file() for p in host_paths):
        raise ValueError("Incomplete or indirect host process evidence")
    evidence.extend(record(p) for p in host_paths)
    summary = json.loads(host_paths[1].read_text())
    if not same(summary, measured["host_process_observation"]):
        raise ValueError("Embedded host process summary differs")
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
    job_paths = [directory / f"job_memory_{name}.json" for name in ("before", "after")]
    if any(p.is_symlink() or not p.is_file() for p in job_paths):
        raise ValueError("Require complete job memory evidence")
    evidence.extend(record(p) for p in job_paths)
    before, after = [json.loads(p.read_text()) for p in job_paths]
    if not same(measured["job_memory"], dict(before=before, after=after)):
        raise ValueError("Embedded job memory differs")
    job_memory = validate_job_memory(before, after, job_scope, done, memory)
    path = directory / "report_finalization.json"
    if path.is_symlink() or not path.is_file():
        raise ValueError("Missing direct reporting-stage evidence")
    evidence.append(record(path))
    finalization = json.loads(path.read_text())
    final_memory = validate_finalization(finalization, job_id, after)
    validate_job_memory(before, final_memory, job_scope, done, memory)
    old, new = (counters(row["raw"]["memory.events"]) for row in (after, final_memory))
    if (int(final_memory["raw"]["memory.peak"]) < int(after["raw"]["memory.peak"])
            or old.keys() != new.keys() or any(new[k] < v for k, v in old.items())):
        raise ValueError("Reporting-stage memory peak or events decreased")
    reproduced = dict(status="native_root_context_measured", job_id=job_id, native_wall_s=wall,
        context=evaluate(points, job_id), lineage_report=lineage_identity(directory),
        scientific_timings_admitted=False, environmental_validity_established=False)
    if not same(context, reproduced):
        raise ValueError("Supplementary context does not reproduce")
    for item in evidence:
        check(item)
    if entries != sorted(directory.iterdir()) or point_files != point_inventory(directory):
        raise ValueError("Raw point inventory changed during replay")
    return dict(schema="allocated_threadripper_replay_v1",
        status="threadripper_scaling_measurement_replayed", native_outcome=outcome,
        native_exit_code=done["exit_code"], native_wall_s=wall, measured=measured, memory=memory,
        host_process_replay=host_replay, affinity_observation_statuses=affinity,
        native_cpu_ids=allowed, allocated_placement=ready["allocated_placement"],
        job_memory=job_memory, job_memory_required=True, report_finalization=finalization,
        native_completion=completion, memory_events=events, screening=screening,
        context=reproduced["context"],
        original_flagged_intervals=screening["original_screening"]["original_threshold_screen"]["flagged_intervals"],
        narrow_flagged_intervals=screening["narrow_flagged_intervals"], evidence=evidence,
        source=record(__file__), common_replay_source=record(Path(__file__).with_name("replay_threadripper_scaling.py")),
        scientific_timings_admitted=False, environmental_validity_established=False,
        native_outputs_validated=False, publication_ready=False,
        limitations=["Replay reproduces allocated placement and accounting, not scientific output validity.",
            "Separate scheduler/session, runtime, frozen recipe and environment audits remain mandatory.",
            "Periodic affinity/host evidence is not continuous compliance or host isolation.",
            "Preserve all outcomes, raw windows and flags; no overhead subtraction or fastest-run selection."])

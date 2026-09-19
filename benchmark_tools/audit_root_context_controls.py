"""Bind terminal scheduler, recipe, runtime and every root-context control trial."""

import argparse
import json
import math
from pathlib import Path
import re

from benchmark_tools.audit_frontier_overhead import recipe_evidence
from benchmark_tools.capture_job_scheduler import terminal_record
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.replay_full_node_control import inventory, read
from benchmark_tools.replay_root_context_control import replay
from benchmark_tools.root_context_control_design import unit_name
from benchmark_tools.run_root_context_controls import PROTOCOL_SHA, PLAN_SHA, summarize
from benchmark_tools.verify_dual_native_provenance import same

ROOT = Path("/home/jlsteenwyk/projects/orthohmm-publication")
RECIPE = "root_context_controls_recipe_v1"
OUTPUT = "root_context_controls_v1"
DEPLOYMENTS = {
    "v1": (RECIPE, OUTPUT, "run_dgx_root_context_controls.sh"),
    "v2": ("root_context_controls_recipe_v2", "root_context_controls_v2", "run_dgx_root_context_session.sh"),
}


def scheduler(raw, job, deployment="v1"):
    recipe_name, _, script = DEPLOYMENTS[deployment]
    unit_name(job, 0)
    retained = terminal_record(raw, job)
    if retained is None:
        raise ValueError("Require terminal root-control scheduler record")
    fields = dict(re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", retained))
    expected = dict(JobId=str(job), Restarts="0", Requeue="0", NodeList="spark-7ff0", Partition="spark",
        OverSubscribe="NO", MinMemoryNode="96G", NumNodes="1", NumCPUs="20", NumTasks="1",
        TimeLimit="00:15:00", WorkDir=str(ROOT / recipe_name),
        Command=str(ROOT / recipe_name / "benchmark_tools" / script))
    expected["CPUs/Task"] = "20"
    if any(fields.get(key) != value for key, value in expected.items()):
        raise ValueError("Root-control scheduler allocation or command differs")
    return fields


def bind(plan, recipe, recipe_sha, verification, launch, panel, job, deployment="v1"):
    recipe_name, _, _ = DEPLOYMENTS[deployment]
    files = {row["path"]: row for row in recipe["records"] if row["kind"] == "file"}
    if len(files) != sum(row["kind"] == "file" for row in recipe["records"]):
        raise ValueError("Duplicate recipe file")
    base = ROOT / recipe_name / "benchmark_tools"
    if files[str(base / "results/ROOT_CPU_CONTEXT_CONTROL_PROTOCOL_20260919.md")]["sha256"] != PROTOCOL_SHA:
        raise ValueError("Frozen protocol differs")
    expected = [dict(row, records=count, status="runtime_tree_identity_matches", scientific_execution_authorized=False)
        for row, count in zip(plan["runtime_manifests"], (26673, 10066))]
    expected.append(dict(path=str(ROOT / (recipe_name + ".json")), sha256=recipe_sha,
        records=len(recipe["records"]), status="runtime_tree_identity_matches", scientific_execution_authorized=False))
    if (len(plan["runtime_manifests"]) != 2 or not same(verification["before"], expected)
            or not same(verification["after"], expected) or not same(verification["measurement"], panel)
            or verification["status"] != panel["status"]
            or verification["source_sha256"] != files[str(base / "run_verified_slurm_measurement.py")]["sha256"]
            or verification["scientific_results_admitted"] is not False):
        raise ValueError("Runtime checks or embedded panel differ")
    for key in ("before_check_wall_s", "after_check_wall_s"):
        value = verification[key]
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            raise ValueError("Invalid runtime-check duration")
    expected_launch = dict(job_id=job, recipe_sha256=recipe_sha, protocol_sha256=PROTOCOL_SHA,
        runtime_plan_sha256=PLAN_SHA, executable=plan["launcher_python"], host="spark-7ff0",
        manager=panel["manager"], source_directory=str(ROOT / recipe_name), scientific_timings_admitted=False)
    if any(not same(launch.get(key), value) for key, value in expected_launch.items()):
        raise ValueError("Launch identity differs")
    if (panel["job_id"] != job or panel["protocol_sha256"] != PROTOCOL_SHA
            or panel["scientific_timings_admitted"] is not False or panel["publication_ready"] is not False
            or not same(panel["summary"], summarize(panel["trials"]))):
        raise ValueError("Panel summary or identity differs")
    failed = None
    for index, row in enumerate(panel["trials"]):
        if row["job_id"] != job or row["scientific_timings_admitted"] is not False:
            raise ValueError("Trial identity differs")
        if failed is not None:
            if row["status"] != "not_run_after_failure" or row.get("failed_index") != failed:
                raise ValueError("Work launched after panel failure")
        elif row["status"] == "failed":
            failed = index
        elif row["status"] != "root_context_workload_validated":
            raise ValueError("Unexpected trial status")
    expected_status = "root_context_controls_completed" if failed is None else "root_context_controls_stopped"
    if panel["status"] != expected_status:
        raise ValueError("Panel terminal status differs")


def audit(archive, recipe_path, recipe_sha, scheduler_path, job, deployment="v1"):
    recipe_name, output_name, _ = DEPLOYMENTS[deployment]
    scheduler_source, recipe_source = record(scheduler_path), record(recipe_path)
    allocation = scheduler(scheduler_path.read_text(), job, deployment)
    recipe = read_pinned(recipe_path, recipe_sha)
    sources = recipe_evidence(archive, recipe, recipe_name)
    plan_path = archive / recipe_name / "benchmark_tools/results/dgx_dual_native_plan_20260919.json"
    plan = read_pinned(plan_path, PLAN_SHA)
    directory = archive / output_name
    evidence = inventory(directory)
    verification, launch = [read(directory / name) for name in ("verification.json", "launch.json")]
    panel = read(directory / "measurement/result.json")
    bind(plan, recipe, recipe_sha, verification, launch, panel, job, deployment)
    expected_exit = "0:0" if (panel["summary"]["all_workloads_valid"]
        and panel["summary"]["all_positive_controls_detected"]) else "1:0"
    if allocation["ExitCode"] != expected_exit or allocation["JobState"] != ("COMPLETED" if expected_exit == "0:0" else "FAILED"):
        raise ValueError("Scheduler result disagrees with panel outcome")
    measurement = directory / "measurement"
    expected_names = {"result.json"} | {f"{prefix}_{i:02d}{suffix}" for i in range(12)
        for prefix, suffix in (("trial", ""), ("progress", ".json"))}
    if {path.name for path in measurement.iterdir()} != expected_names:
        raise ValueError("Complete panel inventory differs")
    rows, failed, previous_end, identity = [], None, None, None
    for index, original in enumerate(panel["trials"]):
        path = measurement / f"trial_{index:02d}"
        if not same(read(path / "panel_trial.json"), original):
            raise ValueError("Retained trial differs from panel")
        if original["status"] == "failed":
            failed = index
        checkpoint = dict(trials=panel["trials"][:index+1], job_id=job, failed_index=failed)
        if not same(read(measurement / f"progress_{index:02d}.json"), checkpoint):
            raise ValueError("Cumulative checkpoint differs")
        row = dict(index=index, retained=original)
        if original["status"] == "root_context_workload_validated":
            result = replay(path, ROOT / output_name / "measurement" / path.name, ROOT / recipe_name,
                            plan["launcher_python"], original["mode"], job, index, panel["manager"])
            if not same(dict(result["trial"], block=original["block"]), original):
                raise ValueError("Replayed trial differs from panel")
            points = result["measurement"]["lineage"]["measured"]["points"]
            first, last = points[0]["root_context"], points[-1]["root_context"]
            current_identity = {key: first[key] for key in ("boot_before", "ticks", "identities_before")}
            if identity is not None and not same(current_identity, identity):
                raise ValueError("Boot or named-scope identity changed across panel")
            identity = current_identity
            start = points[0]["host"][0]["started_monotonic_ns"]
            if previous_end is not None and start <= previous_end:
                raise ValueError("Trial observation or owned-service intervals overlap")
            previous_end = last["host_after"]["finished_ns"]
            if original["mode"] == "user-contended":
                previous_end = max(previous_end, read(path / "controller.json")["removal"][-1]["monotonic_ns"])
            row.update(status="validated", replay=result)
        elif original["status"] == "not_run_after_failure":
            if {p.name for p in path.iterdir()} != {"panel_trial.json"}:
                raise ValueError("Unrun condition contains workload evidence")
            row["status"] = "retained_unrun"
        else:
            row["status"] = "retained_failure"
        rows.append(row)
    for item in [scheduler_source, recipe_source, *sources, *evidence]:
        check(item)
    if inventory(directory) != evidence:
        raise ValueError("Panel evidence changed during audit")
    return dict(status="root_context_controls_audited", trials=rows, allocation=allocation,
        validated_trials=sum(row["status"] == "validated" for row in rows),
        scheduler_source=scheduler_source, recipe_source=recipe_source, recipe_files=sources, evidence=evidence,
        source=record(__file__), scientific_timings_admitted=False, publication_ready=False,
        limitations=["Failed and unrun conditions remain failures or missing observations, not validated controls.",
            "Queue-before-submission evidence and descriptive condition comparisons require separate checks.",
            "No native overhead, process attribution or scientific resource comparison is established."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("archive", "recipe", "scheduler", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--recipe-sha", required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--deployment", choices=DEPLOYMENTS, default="v1")
    args = parser.parse_args()
    result = audit(args.archive.resolve(), args.recipe.resolve(), args.recipe_sha, args.scheduler.resolve(), args.job, args.deployment)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")

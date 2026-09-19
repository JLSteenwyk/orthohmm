"""Audit the terminal full-node control panel without scientific admission."""

import argparse
import json
from pathlib import Path
import re

from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.run_full_node_controls import ORDER, PROTOCOL_SHA, PLAN_SHA, summarize
from benchmark_tools.replay_full_node_control import replay, inventory, read
from benchmark_tools.audit_frontier_overhead import recipe_evidence
from benchmark_tools.capture_job_scheduler import terminal_record
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.verify_dual_native_provenance import same

ROOT = Path("/home/jlsteenwyk/projects/orthohmm-publication")
RECIPE = "full_node_controls_recipe_v1"
RECIPE_SHA = "796da266cf839441d239b1ffae4c372ff616b8feacd64dada84746f0d6c31e62"
JOB = 21918


def scheduler(text):
    if terminal_record(text, JOB) is None:
        raise ValueError("Control panel must be terminal before archive inspection")
    fields = dict(re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", text))
    expected = dict(JobId=str(JOB), Restarts="0", Requeue="0", NodeList="spark-7ff0",
                    OverSubscribe="NO", MinMemoryNode="96G", NumNodes="1", NumCPUs="20",
                    NumTasks="1", TimeLimit="00:45:00", WorkDir=str(ROOT / RECIPE))
    expected["CPUs/Task"] = "20"
    if any(fields.get(key) != value for key, value in expected.items()):
        raise ValueError("Control scheduler allocation differs")
    return fields


def bind(plan, recipe, verification, launch, panel):
    expected = [dict(row, records=count, status="runtime_tree_identity_matches",
                     scientific_execution_authorized=False)
                for row, count in zip(plan["runtime_manifests"], (26673, 10066))]
    expected.append(dict(path=str(ROOT / (RECIPE + ".json")), sha256=RECIPE_SHA,
        records=len(recipe["records"]), status="runtime_tree_identity_matches", scientific_execution_authorized=False))
    wrapper = next(row["sha256"] for row in recipe["records"]
                   if row["path"] == str(ROOT / RECIPE / "benchmark_tools/run_verified_slurm_measurement.py"))
    if (len(plan["runtime_manifests"]) != 2 or not same(verification["before"], expected)
            or not same(verification["after"], expected) or not same(verification["measurement"], panel)
            or verification["status"] != "full_node_controls_completed"
            or verification["source_sha256"] != wrapper or verification["scientific_results_admitted"] is not False):
        raise ValueError("Runtime verification or embedded panel differs")
    expected_launch = dict(job_id=JOB, recipe_sha256=RECIPE_SHA, protocol_sha256=PROTOCOL_SHA,
                           plan_sha256=PLAN_SHA, executable=plan["launcher_python"], host="spark-7ff0")
    if any(not same(launch[key], value) for key, value in expected_launch.items()):
        raise ValueError("Launch identity differs")
    if (panel["job_id"] != JOB or panel["protocol_sha256"] != PROTOCOL_SHA
            or panel["scientific_timings_admitted"] is not False or panel["publication_ready"] is not False
            or not same(panel["summary"], summarize(panel["trials"]))):
        raise ValueError("Panel summary or provenance differs")


def audit(archive, results, scheduler_path):
    scheduler_source = record(scheduler_path)
    allocation = scheduler(scheduler_path.read_text())
    plan = read_pinned(results / "dgx_dual_native_plan_20260919.json", PLAN_SHA)
    recipe = read_pinned(results / "full_node_controls_recipe_20260919.json", RECIPE_SHA)
    recipe_files = recipe_evidence(archive, recipe, RECIPE)
    directory = archive / "full_node_controls_v1"
    evidence = inventory(directory)
    verification, launch = [read(directory / name) for name in ("verification.json", "launch.json")]
    panel = read(directory / "measurement/result.json")
    bind(plan, recipe, verification, launch, panel)
    rows = []
    for index, original in enumerate(panel["trials"]):
        path = directory / "measurement" / f"trial_{index:02d}"
        if not same(original, read(path / "panel_trial.json")):
            raise ValueError("Panel trial differs from retained row")
        row = dict(index=index, block=original["block"], mode=original["mode"], retained=original)
        if original["status"] == "workload_validated":
            try:
                result = replay(path, ROOT / "full_node_controls_v1/measurement" / path.name,
                                ROOT / RECIPE, plan["launcher_python"], original["mode"], JOB)
                if not same(dict(result["trial"], block=original["block"]), original):
                    raise ValueError("Replayed trial differs from panel")
                measured = result["measurement"]
                row.update(status="validated", replay=result,
                    original_flags=len(measured["original_flagged_intervals"]),
                    narrow_flags=len(measured["narrow_flagged_intervals"]),
                    intervals=len(measured["screening"]["narrow_intervals"]),
                    common_intervals=len(result["trial"]["common_intervals"]),
                    common_flags=len(result["trial"]["common_narrow_flagged"]),
                    positive_control_detected=result["trial"]["positive_control_detected"])
            except (OSError, ValueError, KeyError, TypeError) as error:
                row.update(status="invalid", error_type=type(error).__name__, error=str(error))
        else:
            row["status"] = "retained_failure"
        rows.append(row)
    for item in [scheduler_source, *recipe_files, *evidence]:
        check(item)
    if inventory(directory) != evidence:
        raise ValueError("Panel evidence changed during audit")
    return dict(status="full_node_controls_audited", trials=rows,
        validated_trials=sum(row["status"] == "validated" for row in rows),
        allocation=allocation, scheduler_source=scheduler_source, recipe_files=recipe_files,
        evidence=evidence, source=record(__file__), scientific_timings_admitted=False,
        publication_ready=False, limitations=["No threshold changes, timing corrections or selective retries.",
        "Control validity and detection are distinct; neither identifies historical foreign activity.",
        "Scheduler capture replay and descriptive within-block analysis remain separate checks."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("archive", "results", "scheduler", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.archive.resolve(), args.results.resolve(), args.scheduler.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")

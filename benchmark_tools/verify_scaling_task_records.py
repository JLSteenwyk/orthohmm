"""Bind prepared scaling work and recorded runtime checks, not environmental admission."""

from copy import deepcopy
import json
import math
from pathlib import Path

from benchmark_tools.gnu_time_companion import command as time_command
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.measure_root_context_scaling import load_task, PLAN_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_root_context_scaling import ROOT
from benchmark_tools.replay_scaling_root_context import native_outcome
from benchmark_tools.verify_lineage_native_provenance import same

RECIPE_ROOT = ROOT / "root_context_scaling_recipe_v1"
RECIPE_PATH = ROOT / "root_context_scaling_recipe_v1.json"


def bind(plan_path, index, directory, recipe_path, recipe_sha, job):
    if type(job) is not int or job <= 0:
        raise ValueError("Require explicit positive job identity")
    directory = Path(directory)
    paths = [directory / name for name in (
        "preparation.json", "verification.json", "measurement/lineage_report.json", "measurement/root_context_report.json")]
    if directory.is_symlink() or any(p.is_symlink() for p in paths):
        raise ValueError("Direct task evidence links not supported")
    sources = [record(plan_path), record(recipe_path)]
    if sources[0]["sha256"] != PLAN_SHA or sources[1]["sha256"] != recipe_sha:
        raise ValueError("Plan or recipe hash differs")
    plan, task, order = load_task(plan_path, index)
    recipe = read_pinned(recipe_path, recipe_sha)
    entries = recipe["records"]
    if (recipe["roots"] != [str(RECIPE_ROOT)]
            or len({r["path"] for r in entries}) != len(entries)):
        raise ValueError("Unexpected recipe root or duplicate records")
    file_paths = {r["path"] for r in entries if r["kind"] == "file"}
    if any(str(RECIPE_ROOT / "benchmark_tools" / name) not in file_paths for name in (
            "run_verified_slurm_measurement.py", "measure_scaling_root_context.py",
            "measure_root_context_scaling.py", "replay_scaling_root_context.py")):
        raise ValueError("Missing required measurement source in recipe")
    evidence = [record(p) for p in paths]
    preparation, verification, measured, context = [json.loads(p.read_text()) for p in paths]
    if type(measured["job_id"]) is not int or measured["job_id"] != job:
        raise ValueError("Measurement job identity differs")
    outcome, wall = native_outcome(measured, measured["native"])
    if not measured["points"] or any("root_context" not in p for p in measured["points"]):
        raise ValueError("Missing supplementary observation fields")
    if not same(context["lineage_report"], dict(evidence[2], path="lineage_report.json")):
        raise ValueError("Supplementary report binds different lineage bytes")
    run = deepcopy(task["run"])
    remote = Path(run["measurement_directory"]).parent
    run["gnu_time"] = dict(executable="/usr/bin/time", output=str(remote / "native.time.tsv"))
    argv = time_command(run["native_argv"], run["gnu_time"]["output"])
    originals, copies = order["inputs_in_native_order"], []
    if run["native_method"] == "orthofinder_full":
        copies = [dict(row, path=str(Path(run["configuration"]["copy_inputs_to"]) / Path(row["path"]).name))
                  for row in sorted(originals, key=lambda r:Path(r["path"]).name)]
    expected = dict(run=run, measured_argv=argv, original_inputs=originals, copied_inputs=copies,
        expected_native_basename_order=sorted(order["native_order"]) if copies else order["native_order"],
        status="fresh_native_inputs_prepared", inference_started=False)
    if any(not same(preparation[k], v) for k,v in expected.items()):
        raise ValueError("Prepared work, input bytes or enumeration differ")
    if len(plan["runtime_manifests"]) != 2:
        raise ValueError("Unexpected runtime inventory")
    runtime = [dict(row, records=count, scientific_execution_authorized=False, status="runtime_tree_identity_matches")
               for row,count in zip(plan["runtime_manifests"], (26673,10066))]
    runtime.append(dict(path=str(RECIPE_PATH), sha256=recipe_sha, records=len(entries),
                        scientific_execution_authorized=False, status="runtime_tree_identity_matches"))
    before = dict(runtime=runtime, original_inputs=originals, native_order=order["native_order"])
    after = deepcopy(before)
    if copies:
        after["copied_inputs"] = copies
    wrappers = [r["sha256"] for r in entries if r["kind"] == "file"
                and r["path"] == str(RECIPE_ROOT / "benchmark_tools/run_verified_slurm_measurement.py")]
    if len(wrappers) != 1:
        raise ValueError("Missing or ambiguous wrapper source")
    expected = dict(status=measured["status"], before=before, after=after, measurement=measured,
                    source_sha256=wrappers[0], scientific_results_admitted=False)
    if any(not same(verification[k], v) for k,v in expected.items()):
        raise ValueError("Recorded runtime, wrapper or embedded measurement differs")
    for value in (preparation["preparation_wall_s"], verification["before_check_wall_s"], verification["after_check_wall_s"]):
        if type(value) not in (int,float) or not math.isfinite(value) or value <= 0:
            raise ValueError("Invalid preparation/verification duration")
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
        plan["launcher_python"], "-B", str(RECIPE_ROOT / "benchmark_tools/measure_scaling_root_context.py"),
        "--worker", run["measurement_directory"]]
    if not same(measured["launched"], launched):
        raise ValueError("Worker entry point or native step command differs")
    for item in [*sources, *evidence]:
        check(item)
    return dict(status="scaling_task_records_bound", index=index, job_id=job, task=task, run=run,
        measured_argv=argv, measurement=measured, verification=verification, native_outcome=outcome,
        native_wall_s=wall, sources=sources, evidence=evidence, source=record(__file__),
        scientific_timings_admitted=False, native_outputs_validated=False, environmental_validity_established=False,
        limitations=["Recorded task/runtime claims bound, not independent current runtime-tree verification.",
            "Requires recipe archive, authorization, scheduler/session, raw replay and native-output/failure audits.",
            "Before/after identities cannot exclude temporary changes during execution; no timing admission."])

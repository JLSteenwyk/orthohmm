"""Bind a successful native root-context task; raw/output audit remains separate."""

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import re

from benchmark_tools.gnu_time_companion import command as time_command
from benchmark_tools.capture_job_scheduler import terminal_record
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_root_context_native import ROOT, build
from benchmark_tools.run_root_context_native import PLAN_SHA
from benchmark_tools.verify_lineage_native_provenance import same

RECIPE_ROOT = ROOT / "root_context_native_recipe_v1"
RECIPE_PATH = ROOT / "root_context_native_recipe_v1.json"
SUBMISSION_SCRIPT = RECIPE_ROOT / "benchmark_tools/run_dgx_root_context_native.sh"


def load_context(results, recipe_path, recipe_sha):
    if not re.fullmatch(r"[0-9a-f]{64}", recipe_sha):
        raise ValueError("Invalid recipe hash")
    plan = read_pinned(results / "dgx_root_context_native_plan_20260919.json", PLAN_SHA)
    expected = build(results / "dgx_lineage_native_plan_20260919.json",
                     results / "ROOT_CONTEXT_NATIVE_PROTOCOL_20260919.md")
    if not same(plan, expected):
        raise ValueError("Frozen plan derivation differs")
    recipe = read_pinned(recipe_path, recipe_sha)
    paths = [r["path"] for r in recipe["records"]]
    if len(paths) != len(set(paths)):
        raise ValueError("Duplicate recipe records")
    return dict(plan=plan, recipe=recipe, recipe_sha256=recipe_sha)


def scheduler_identity(text, job, require_completed=True):
    if type(job) is not int or job <= 0:
        raise ValueError("Invalid job identity")
    fields = {}
    for key, value in re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", text):
        if key in fields:
            raise ValueError("Duplicate scheduler field")
        fields[key] = value
    expected = dict(JobId=str(job), JobState="COMPLETED", ExitCode="0:0", Restarts="0",
        Requeue="0", NodeList="spark-7ff0", OverSubscribe="NO", MinMemoryNode="96G",
        NumNodes="1", NumCPUs="20", NumTasks="1", TimeLimit="01:00:00",
        Command=str(SUBMISSION_SCRIPT), WorkDir=str(RECIPE_ROOT))
    expected["CPUs/Task"] = "20"
    if not require_completed:
        if terminal_record(text, job) is None:
            raise ValueError("Require terminal allocation")
        del expected["JobState"]
        del expected["ExitCode"]
    if any(fields.get(k) != v for k, v in expected.items()) or any(
            k in fields for k in ("ArrayJobId", "ArrayTaskId")):
        raise ValueError("Scheduler identity, resources or completion differs")
    return job


def verify(context, index, preparation, verification_path, receipt, measurement, scheduler, job,
           require_completed=True):
    scheduler_identity(scheduler, job, require_completed)
    if type(index) is not int or index not in range(3):
        raise ValueError("Invalid task index")
    plan, recipe = context["plan"], context["recipe"]
    task = plan["runs"][index]
    raw = Path(verification_path).read_bytes()
    verification = json.loads(raw)
    run = deepcopy(task["run"])
    directory = Path(run["measurement_directory"])
    expected_receipt = dict(index=index, task=task, job_id=job,
        scientific_timings_admitted=False, status="measurement_completed",
        wrapper_status="command_exited_zero", native_wall_s=measurement["native_wall_s"],
        verification=dict(path=str(directory.parent / "verification.json"),
                          bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    if not same(receipt, expected_receipt):
        raise ValueError("Task receipt or verification bytes differ")
    if type(measurement["job_id"]) is not int or measurement["job_id"] != job:
        raise ValueError("Measurement job differs")
    native = measurement["native"]
    if type(native["exit_code"]) is not int or native["exit_code"] != 0 or native["timed_out"] is not False:
        raise ValueError("Native command did not succeed")
    run["gnu_time"] = dict(executable="/usr/bin/time", output=str(directory.parent / "native.time.tsv"))
    argv = time_command(run["native_argv"], run["gnu_time"]["output"])
    order = plan["order"]
    originals = order["inputs_in_native_order"]
    copies = []
    if run["native_method"] == "orthofinder_full":
        copies = [dict(item, path=str(Path(run["configuration"]["copy_inputs_to"]) / Path(item["path"]).name))
                  for item in sorted(originals, key=lambda item: Path(item["path"]).name)]
    expected = dict(run=run, measured_argv=argv, original_inputs=originals, copied_inputs=copies,
        expected_native_basename_order=sorted(order["native_order"]) if copies else order["native_order"],
        status="fresh_native_inputs_prepared", inference_started=False)
    if any(not same(preparation[k], v) for k, v in expected.items()):
        raise ValueError("Preparation differs from frozen task")
    if len(plan["runtime_manifests"]) != 2:
        raise ValueError("Unexpected runtime inventory")
    runtime = [dict(item, records=count, scientific_execution_authorized=False,
                    status="runtime_tree_identity_matches")
               for item, count in zip(plan["runtime_manifests"], (26673, 10066))]
    runtime.append(dict(path=str(RECIPE_PATH), sha256=context["recipe_sha256"],
        records=len(recipe["records"]), scientific_execution_authorized=False,
        status="runtime_tree_identity_matches"))
    before = dict(native_order=order["native_order"], original_inputs=originals, runtime=runtime)
    after = deepcopy(before)
    if copies:
        after["copied_inputs"] = copies
    wrappers = [r["sha256"] for r in recipe["records"] if r["kind"] == "file" and
                r["path"] == str(RECIPE_ROOT / "benchmark_tools/run_verified_slurm_measurement.py")]
    if len(wrappers) != 1:
        raise ValueError("Missing or duplicate wrapper")
    expected = dict(status="command_exited_zero", before=before, after=after,
        measurement=measurement, source_sha256=wrappers[0], scientific_results_admitted=False)
    if any(not same(verification[k], v) for k, v in expected.items()):
        raise ValueError("Runtime, wrapper or embedded measurement differs")
    for value in (preparation["preparation_wall_s"], verification["before_check_wall_s"],
                  verification["after_check_wall_s"], measurement["native_wall_s"]):
        if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
            raise ValueError("Invalid duration")
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
        plan["launcher_python"], "-B", str(RECIPE_ROOT / "benchmark_tools/measure_native_lineage_step.py"),
        "--worker", run["measurement_directory"]]
    if not same(measurement["launched"], launched):
        raise ValueError("Worker command differs")
    return dict(status="root_context_native_provenance_bound", job_id=job, index=index,
        run=run, measured_argv=argv, scientific_timings_admitted=False,
        limitations=["Requires separately pinned context, whole-panel and session receipt validation.",
                     "Raw collector replay and native output equivalence are separate checks.",
                     "Before/after runtime verification cannot exclude temporary changes."])

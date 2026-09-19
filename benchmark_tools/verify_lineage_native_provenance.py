"""Bind lineage-native diagnostics to pinned tasks, inputs, runtime and scheduler."""

from copy import deepcopy
import json
import math
from pathlib import Path
import re

from benchmark_tools.gnu_time_companion import command as time_command
from benchmark_tools.run_lineage_native_diagnostic import ROOT, build, read_pinned

PLAN_SHA = "714ef04458904e1c01526d171d0b2fde658ac135bd416af49ab107dc5ed14bfe"
RECIPE_SHA = "d05c1028f486d5e0a9788912fc8d3785a7a2ca504e9caafbe5e8094a7a7dc24b"
JOBS = (21995, 21996, 21997)
RECIPE_ROOT = ROOT / "lineage_native_recipe_v1"
SUBMISSION_SCRIPT = Path("/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmark_tools/run_dgx_lineage_native.sh")


def same(left, right):
    return json.dumps(left, sort_keys=True, allow_nan=False) == json.dumps(right, sort_keys=True, allow_nan=False)


def load_context(results):
    plan = read_pinned(results / "dgx_lineage_native_plan_20260919.json", PLAN_SHA)
    expected = build(results / "dgx_pressure_overhead_plan_v2_20260919.json",
                     results / "LINEAGE_NATIVE_PROTOCOL_20260919.md")
    if not same(plan, expected):
        raise ValueError("Frozen plan differs from prescribed derivation")
    return dict(plan=plan, recipe=read_pinned(results / "dgx_lineage_native_recipe_20260919.json", RECIPE_SHA))


def scheduler_identity(text, index):
    if type(index) is not int or index not in range(3):
        raise ValueError("Invalid diagnostic index")
    fields = {}
    for key, value in re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", text):
        if key in fields:
            raise ValueError("Duplicate scheduler field: " + key)
        fields[key] = value
    expected = dict(JobId=str(JOBS[index]), JobState="COMPLETED", ExitCode="0:0",
                    Restarts="0", Requeue="0", NodeList="spark-7ff0", OverSubscribe="NO",
                    MinMemoryNode="96G", NumNodes="1", NumCPUs="20", NumTasks="1",
                    Command=str(SUBMISSION_SCRIPT))
    expected["CPUs/Task"] = "20"
    if any(fields.get(key) != value for key, value in expected.items()) or any(
            key in fields for key in ("ArrayJobId", "ArrayTaskId")):
        raise ValueError("Scheduler identity, allocation, command or completion differs")
    return JOBS[index]


def verify(context, index, preparation, verification, receipt, measurement, scheduler):
    job = scheduler_identity(scheduler, index)
    plan, recipe = context["plan"], context["recipe"]
    task = plan["runs"][index]
    expected = dict(task=task, plan_sha256=PLAN_SHA, recipe_sha256=RECIPE_SHA,
                    purpose="three_native_lineage_diagnostics", scientific_timings_admitted=False)
    if not same(receipt, expected):
        raise ValueError("Diagnostic task receipt differs")
    if type(measurement["job_id"]) is not int or measurement["job_id"] != job:
        raise ValueError("Measurement and scheduler identities differ")
    run = deepcopy(task["run"])
    run["gnu_time"] = dict(executable="/usr/bin/time",
                          output=str(Path(run["measurement_directory"]).parent / "native.time.tsv"))
    argv = time_command(run["native_argv"], run["gnu_time"]["output"])
    order = plan["order"]
    originals = order["inputs_in_native_order"]
    copied = []
    if run["native_method"] == "orthofinder_full":
        copied = [dict(item, path=str(Path(run["configuration"]["copy_inputs_to"]) / Path(item["path"]).name))
                  for item in sorted(originals, key=lambda item: Path(item["path"]).name)]
    expected_preparation = dict(run=run, measured_argv=argv, original_inputs=originals,
        copied_inputs=copied, expected_native_basename_order=sorted(order["native_order"]) if copied else order["native_order"],
        status="fresh_native_inputs_prepared", inference_started=False)
    if any(not same(preparation[key], value) for key, value in expected_preparation.items()):
        raise ValueError("Prepared command, inputs or order differ from plan")
    if len(plan["runtime_manifests"]) != 2:
        raise ValueError("Unexpected runtime manifest inventory")
    runtime = [dict(item, records=count, scientific_execution_authorized=False, status="runtime_tree_identity_matches")
               for item, count in zip(plan["runtime_manifests"], (26673, 10066))]
    runtime.append(dict(path=str(ROOT / "lineage_native_recipe_v1.json"), sha256=RECIPE_SHA,
                        records=len(recipe["records"]), scientific_execution_authorized=False,
                        status="runtime_tree_identity_matches"))
    before = dict(native_order=order["native_order"], original_inputs=originals, runtime=runtime)
    after = dict(before)
    if copied:
        after["copied_inputs"] = copied
    wrappers = [r["sha256"] for r in recipe["records"]
                if r["path"] == str(RECIPE_ROOT / "benchmark_tools/run_verified_slurm_measurement.py") and r["kind"] == "file"]
    if len(wrappers) != 1:
        raise ValueError("Missing or duplicate frozen wrapper")
    if (verification["status"] != "command_exited_zero" or not same(verification["before"], before)
            or not same(verification["after"], after) or not same(verification["measurement"], measurement)
            or verification["source_sha256"] != wrappers[0] or verification["scientific_results_admitted"] is not False):
        raise ValueError("Runtime, wrapper or embedded measurement verification differs")
    for value in (preparation["preparation_wall_s"], verification["before_check_wall_s"], verification["after_check_wall_s"]):
        if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
            raise ValueError("Invalid preparation or verification duration")
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
                plan["launcher_python"], "-B", str(RECIPE_ROOT / "benchmark_tools/measure_native_lineage_step.py"),
                "--worker", run["measurement_directory"]]
    if not same(measurement["launched"], launched):
        raise ValueError("Wrong native worker/collector launch")
    return dict(status="lineage_native_provenance_bound", job_id=job, run=run, measured_argv=argv,
        scientific_timings_admitted=False,
        limitations=["Record binding only; raw replay, archive hashes and native outputs require separate checks.",
                     "Retained before/after runtime checks cannot exclude temporary changes during execution."])

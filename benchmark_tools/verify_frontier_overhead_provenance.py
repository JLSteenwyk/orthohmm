"""Bind successful overhead task records to frozen commands and provenance."""

from copy import deepcopy
import json
import math
from pathlib import Path
import re

from benchmark_tools.gnu_time_companion import command as time_command
from benchmark_tools.run_dgx_frontier_overhead import PLAN_SHA, ROOT, read_pinned

RECIPE_SHA = "64a05b5201e78a9d8d46f302879a49bd5cc470bf7813eb5d16e4c88d79c51edc"
AUTH_SHA = "77875e0273883454c25fcb697916aedc10f5d0d3081baa77ead997d6f4aa0b47"


def load_context(results):
    return {
        "plan": read_pinned(results / "dgx_frontier_overhead_plan_20260918.json", PLAN_SHA),
        "recipe": read_pinned(results / "dgx_frontier_overhead_recipe_v1_20260918.json", RECIPE_SHA),
        "authorization": read_pinned(results / "dgx_frontier_overhead_authorization_20260918.json", AUTH_SHA),
    }


def scheduler_identity(text, index, array_id=21838):
    fields = {}
    for key, value in re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", text):
        if key in fields:
            raise ValueError("Duplicate scheduler field: " + key)
        fields[key] = value
    expected = dict(ArrayJobId=str(array_id), ArrayTaskId=str(index), JobState="COMPLETED",
                    ExitCode="0:0", Restarts="0", Requeue="0", NodeList="spark-7ff0",
                    OverSubscribe="NO", MinMemoryNode="96G", NumNodes="1", NumCPUs="20")
    expected["CPUs/Task"] = "20"
    expected["Command"] = str(ROOT / "frontier_overhead_recipe_v1/benchmark_tools/run_dgx_frontier_overhead.sh")
    if any(fields.get(key) != value for key, value in expected.items()):
        raise ValueError("Scheduler task, allocation, command or completion differs")
    if not fields.get("JobId", "").isdigit() or int(fields["JobId"]) <= 0:
        raise ValueError("Invalid native scheduler job identity")
    return int(fields["JobId"])


def verify(context, index, preparation, verification, receipt, measurement, scheduler):
    if type(index) is not int or not 0 <= index < 18:
        raise ValueError("Invalid task index")
    plan, recipe, authorization = (context[key] for key in ("plan", "recipe", "authorization"))
    task = plan["runs"][index]
    expected_receipt = dict(task=task, authorization=authorization, authorization_sha256=AUTH_SHA,
                            plan_sha256=PLAN_SHA, scientific_timings_admitted=False)
    if json.dumps(receipt, sort_keys=True, allow_nan=False) != json.dumps(expected_receipt, sort_keys=True, allow_nan=False):
        raise ValueError("Task authorization receipt differs")
    job_id = scheduler_identity(scheduler, index)
    if measurement["job_id"] != job_id:
        raise ValueError("Measurement and scheduler job identities differ")
    run = deepcopy(task["run"])
    directory = Path(run["measurement_directory"]).parent
    run["gnu_time"] = dict(executable="/usr/bin/time", output=str(directory / "native.time.tsv"))
    argv = time_command(run["native_argv"], run["gnu_time"]["output"])
    order = plan["order"]
    originals = order["inputs_in_native_order"]
    copied = []
    if run["native_method"] == "orthofinder_full":
        copied = [dict(item, path=str(Path(run["configuration"]["copy_inputs_to"]) / Path(item["path"]).name))
                  for item in sorted(originals, key=lambda item: Path(item["path"]).name)]
    expected_order = sorted(order["native_order"]) if copied else order["native_order"]
    if (preparation["run"] != run or preparation["measured_argv"] != argv
            or preparation["original_inputs"] != originals or preparation["copied_inputs"] != copied
            or preparation["expected_native_basename_order"] != expected_order
            or preparation["status"] != "fresh_native_inputs_prepared" or preparation["inference_started"] is not False):
        raise ValueError("Native command or input preparation differs from plan")
    recipe_root = ROOT / "frontier_overhead_recipe_v1"
    runtime = [dict(item, records=count, scientific_execution_authorized=False, status="runtime_tree_identity_matches")
               for item, count in zip(plan["runtime_manifests"], (26673, 10066))]
    runtime.append(dict(path=str(ROOT / "frontier_overhead_recipe_v1.json"), sha256=RECIPE_SHA,
                        records=len(recipe["records"]), scientific_execution_authorized=False,
                        status="runtime_tree_identity_matches"))
    before = dict(native_order=order["native_order"], original_inputs=originals, runtime=runtime)
    after = dict(before)
    if copied:
        after["copied_inputs"] = copied
    wrapper_path = str(recipe_root / "benchmark_tools/run_verified_slurm_measurement.py")
    wrapper_sha = next(row["sha256"] for row in recipe["records"] if row["path"] == wrapper_path)
    if (verification["status"] != "command_exited_zero" or verification["before"] != before
            or verification["after"] != after or verification["measurement"] != measurement
            or verification["source_sha256"] != wrapper_sha or verification["scientific_results_admitted"] is not False):
        raise ValueError("Runtime, wrapper source or embedded measurement verification differs")
    for value in [preparation["preparation_wall_s"], verification["before_check_wall_s"], verification["after_check_wall_s"]]:
        if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
            raise ValueError("Invalid preparation or verification duration")
    collector = "measure_frontier_boundary_step.py" if task["mode"] == "boundary" else "measure_native_frontier_step.py"
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
                str(ROOT / "envs/orthohmm/bin/python"), "-B", str(recipe_root / "benchmark_tools" / collector),
                "--worker", run["measurement_directory"]]
    if measurement["launched"] != launched:
        raise ValueError("Wrong native worker/collector launch")
    return dict(status="overhead_task_provenance_bound", job_id=job_id, run=run, measured_argv=argv,
                scientific_timings_admitted=False,
                limitations=["Record binding only; raw measurement replay and native-output validation remain separate.",
                             "Retained runtime checks cannot exclude temporary changes during execution."])

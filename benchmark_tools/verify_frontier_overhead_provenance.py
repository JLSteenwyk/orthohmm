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

PANELS = {
    "frontier_21838": dict(array_id=21838, plan_sha=PLAN_SHA, recipe_sha=RECIPE_SHA, auth_sha=AUTH_SHA,
        plan_file="dgx_frontier_overhead_plan_20260918.json",
        recipe_file="dgx_frontier_overhead_recipe_v1_20260918.json",
        auth_file="dgx_frontier_overhead_authorization_20260918.json",
        recipe_root="frontier_overhead_recipe_v1", recipe_manifest="frontier_overhead_recipe_v1.json",
        scheduler_command=str(ROOT / "frontier_overhead_recipe_v1/benchmark_tools/run_dgx_frontier_overhead.sh"),
        protocols=["DGX_NATIVE_FRONTIER_OVERHEAD_PROTOCOL_20260918.md"], pressure=False),
    "pressure_21889": dict(array_id=21889,
        plan_sha="b644e165dbf4d0beabf1cf4d9b6c314de522e3ebd1b91598ebebea99094c8fff",
        recipe_sha="50fc3a14c4d5b53c4eb158a22efcfb8b289b346718fad3c78b300172e319d8a1",
        auth_sha="17943df08dfe21bf1975ff0dcb29dc66a9507b156e4698b1672e7517cdfd85e6",
        plan_file="dgx_pressure_overhead_plan_v2_20260919.json",
        recipe_file="dgx_pressure_overhead_recipe_v2_20260919.json",
        auth_file="dgx_pressure_overhead_authorization_v2_20260919.json",
        recipe_root="pressure_overhead_recipe_v2", recipe_manifest="pressure_overhead_recipe_v2_manifest.json",
        scheduler_command="/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmark_tools/run_dgx_pressure_overhead_v2.sbatch",
        protocols=["DGX_PRESSURE_OVERHEAD_PROTOCOL_20260919.md", "DGX_PRESSURE_OVERHEAD_V2_20260919.md"],
        pressure=True),
    "dual_21920": dict(array_id=21920,
        plan_sha="1160a669a4033c66e7bdde5baddf429e83b9328d3821f99291fd29210cac8ab9",
        recipe_sha="ba18b0dca91698aa07cb451085c9340346828eaf86be90b7de53cfc2badceeb8",
        auth_sha="1179b6fc403b146d43aafb2363379c65c1a176afe442914046ec49d84ace4144",
        plan_file="dgx_dual_overhead_plan_20260919.json",
        recipe_file="dgx_dual_overhead_recipe_20260919.json",
        auth_file="dgx_dual_overhead_authorization_20260919.json",
        recipe_root="dual_overhead_recipe_v1", recipe_manifest="dual_overhead_recipe_v1.json",
        scheduler_command="/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/dual_overhead_deployment_6e222a4/dual_overhead_recipe_v1/benchmark_tools/run_dgx_dual_overhead.sh",
        protocols=["DUAL_COLLECTOR_OVERHEAD_PROTOCOL_20260919.md"], pressure=True, dual=True),
}


def panel_spec(panel):
    if panel not in PANELS:
        raise ValueError("Unknown frozen overhead panel")
    return deepcopy(PANELS[panel])


def load_context(results, panel="frontier_21838"):
    spec = panel_spec(panel)
    return {
        "panel": panel,
        "plan": read_pinned(results / spec["plan_file"], spec["plan_sha"]),
        "recipe": read_pinned(results / spec["recipe_file"], spec["recipe_sha"]),
        "authorization": read_pinned(results / spec["auth_file"], spec["auth_sha"]),
    }


def scheduler_identity(text, index, array_id=21838, *, expected_command=None):
    fields = {}
    for key, value in re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", text):
        if key in fields:
            raise ValueError("Duplicate scheduler field: " + key)
        fields[key] = value
    expected = dict(ArrayJobId=str(array_id), ArrayTaskId=str(index), JobState="COMPLETED",
                    ExitCode="0:0", Restarts="0", Requeue="0", NodeList="spark-7ff0",
                    OverSubscribe="NO", MinMemoryNode="96G", NumNodes="1", NumCPUs="20")
    expected["CPUs/Task"] = "20"
    expected["Command"] = expected_command or PANELS["frontier_21838"]["scheduler_command"]
    if any(fields.get(key) != value for key, value in expected.items()):
        raise ValueError("Scheduler task, allocation, command or completion differs")
    if not fields.get("JobId", "").isdigit() or int(fields["JobId"]) <= 0:
        raise ValueError("Invalid native scheduler job identity")
    return int(fields["JobId"])


def verify(context, index, preparation, verification, receipt, measurement, scheduler):
    if type(index) is not int or not 0 <= index < 18:
        raise ValueError("Invalid task index")
    plan, recipe, authorization = (context[key] for key in ("plan", "recipe", "authorization"))
    spec = panel_spec(context.get("panel", "frontier_21838"))
    task = plan["runs"][index]
    expected_receipt = dict(task=task, authorization=authorization, authorization_sha256=spec["auth_sha"],
                            plan_sha256=spec["plan_sha"], scientific_timings_admitted=False)
    if spec.get("dual"):
        expected_receipt["recipe_sha256"] = spec["recipe_sha"]
    if json.dumps(receipt, sort_keys=True, allow_nan=False) != json.dumps(expected_receipt, sort_keys=True, allow_nan=False):
        raise ValueError("Task authorization receipt differs")
    job_id = scheduler_identity(scheduler, index, spec["array_id"], expected_command=spec["scheduler_command"])
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
    recipe_root = ROOT / spec["recipe_root"]
    runtime = [dict(item, records=count, scientific_execution_authorized=False, status="runtime_tree_identity_matches")
               for item, count in zip(plan["runtime_manifests"], (26673, 10066))]
    runtime.append(dict(path=str(ROOT / spec["recipe_manifest"]), sha256=spec["recipe_sha"],
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
    if spec.get("dual") and task["mode"] == "periodic":
        collector = "measure_native_dual_bracket_step.py"
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
                str(ROOT / "envs/orthohmm/bin/python"), "-B", str(recipe_root / "benchmark_tools" / collector),
                "--worker", run["measurement_directory"]]
    if measurement["launched"] != launched:
        raise ValueError("Wrong native worker/collector launch")
    return dict(status="overhead_task_provenance_bound", job_id=job_id, run=run, measured_argv=argv,
                scientific_timings_admitted=False,
                limitations=["Record binding only; raw measurement replay and native-output validation remain separate.",
                             "Retained runtime checks cannot exclude temporary changes during execution."])

from copy import deepcopy
from pathlib import Path

import pytest

from benchmark_tools import verify_dual_native_provenance as module

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def fixture(index):
    context = module.load_context(RESULTS)
    plan, recipe = context["plan"], context["recipe"]
    task = plan["runs"][index]
    run = deepcopy(task["run"])
    run["gnu_time"] = dict(executable="/usr/bin/time", output=str(Path(run["measurement_directory"]).parent / "native.time.tsv"))
    order = plan["order"]
    originals = deepcopy(order["inputs_in_native_order"])
    copies = []
    if index == 2:
        copies = [dict(item, path=str(Path(run["configuration"]["copy_inputs_to"]) / Path(item["path"]).name))
                  for item in sorted(originals, key=lambda item: Path(item["path"]).name)]
    prepared = dict(run=run, measured_argv=module.time_command(run["native_argv"], run["gnu_time"]["output"]),
        original_inputs=originals, copied_inputs=copies,
        expected_native_basename_order=sorted(order["native_order"]) if copies else order["native_order"],
        status="fresh_native_inputs_prepared", inference_started=False, preparation_wall_s=.1)
    runtime = [dict(item, records=count, scientific_execution_authorized=False, status="runtime_tree_identity_matches")
               for item, count in zip(plan["runtime_manifests"], (26673, 10066))]
    runtime.append(dict(path=str(module.ROOT / "dual_native_recipe_v1.json"), sha256=module.RECIPE_SHA,
        records=len(recipe["records"]), scientific_execution_authorized=False, status="runtime_tree_identity_matches"))
    before = dict(runtime=runtime, original_inputs=originals, native_order=order["native_order"])
    after = deepcopy(before)
    if copies:
        after["copied_inputs"] = deepcopy(copies)
    measured = dict(job_id=module.JOBS[index], launched=["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1",
        "--cpus-per-task=20", plan["launcher_python"], "-B",
        str(module.RECIPE_ROOT / "benchmark_tools/measure_native_dual_bracket_step.py"), "--worker", run["measurement_directory"]])
    wrapper = next(r["sha256"] for r in recipe["records"] if r["path"].endswith("/run_verified_slurm_measurement.py"))
    verified = dict(status="command_exited_zero", before=deepcopy(before), after=after,
        measurement=deepcopy(measured), source_sha256=wrapper, scientific_results_admitted=False,
        before_check_wall_s=1., after_check_wall_s=1.)
    receipt = dict(task=deepcopy(task), plan_sha256=module.PLAN_SHA, recipe_sha256=module.RECIPE_SHA,
        purpose="three_native_dual_bracket_diagnostics", scientific_timings_admitted=False)
    scheduler = (f"JobId={module.JOBS[index]} JobState=COMPLETED ExitCode=0:0 Restarts=0 Requeue=0 "
        "NodeList=spark-7ff0 OverSubscribe=NO MinMemoryNode=96G NumNodes=1 NumCPUs=20 NumTasks=1 CPUs/Task=20 "
        f"Command={module.RECIPE_ROOT}/benchmark_tools/run_dgx_dual_native.sh")
    return [context, index, prepared, verified, receipt, measured, scheduler]


@pytest.mark.parametrize("index", range(3))
def test_all_frozen_tasks_bind(index):
    args = fixture(index)
    result = module.verify(*args)
    assert result["job_id"] == module.JOBS[index]
    assert result["run"] == args[2]["run"]
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["receipt", "receipt_bool", "task", "command", "input", "order", "copy",
    "runtime", "runtime_bool", "after", "wrapper", "measurement", "admitted", "duration", "worker",
    "scheduler_job", "scheduler_cpu", "scheduler_restart", "scheduler_exclusive", "scheduler_duplicate", "scheduler_array"])
def test_changed_provenance_rejected(fault):
    args = fixture(2)
    _, _, prepared, verified, receipt, measured, scheduler = args
    if fault == "receipt":
        receipt["recipe_sha256"] = "0"*64
    elif fault == "receipt_bool":
        receipt["scientific_timings_admitted"] = 0
    elif fault == "task":
        receipt["task"]["index"] = 0
    elif fault == "command":
        prepared["run"]["native_argv"].append("--unexpected")
    elif fault == "input":
        prepared["original_inputs"][0]["sha256"] = "0"*64
    elif fault == "order":
        prepared["expected_native_basename_order"] = []
    elif fault == "copy":
        verified["after"].pop("copied_inputs")
    elif fault == "runtime":
        verified["before"]["runtime"][0]["sha256"] = "0"*64
    elif fault == "runtime_bool":
        verified["before"]["runtime"][0]["scientific_execution_authorized"] = 0
    elif fault == "after":
        verified["after"]["native_order"] = []
    elif fault == "wrapper":
        verified["source_sha256"] = "0"*64
    elif fault == "measurement":
        verified["measurement"]["job_id"] += 1
    elif fault == "admitted":
        verified["scientific_results_admitted"] = True
    elif fault == "duration":
        verified["before_check_wall_s"] = float("nan")
    elif fault == "worker":
        measured["launched"][8] = "/wrong/collector.py"
        verified["measurement"] = deepcopy(measured)
    elif fault == "scheduler_duplicate":
        args[-1] += " JobId=21914"
    elif fault == "scheduler_array":
        args[-1] += " ArrayJobId=21914 ArrayTaskId=0"
    else:
        before, after = {"scheduler_job": ("JobId=21914", "JobId=21912"),
            "scheduler_cpu": ("CPUs/Task=20", "CPUs/Task=19"),
            "scheduler_restart": ("Restarts=0", "Restarts=1"),
            "scheduler_exclusive": ("OverSubscribe=NO", "OverSubscribe=YES")}[fault]
        args[-1] = scheduler.replace(before, after)
    with pytest.raises(ValueError):
        module.verify(*args)


@pytest.mark.parametrize("index", [-1, 3, True, 1.0])
def test_invalid_index_rejected(index):
    args = fixture(0)
    args[1] = index
    with pytest.raises(ValueError):
        module.verify(*args)


@pytest.mark.parametrize("name", ["dgx_dual_native_plan_20260919.json", "dgx_dual_native_recipe_20260919.json",
    "dgx_pressure_overhead_plan_v2_20260919.json", "DUAL_BRACKET_NATIVE_PROTOCOL_20260919.md"])
def test_context_pins_source_bytes(tmp_path, name):
    for source in ("dgx_dual_native_plan_20260919.json", "dgx_dual_native_recipe_20260919.json",
                   "dgx_pressure_overhead_plan_v2_20260919.json", "DUAL_BRACKET_NATIVE_PROTOCOL_20260919.md"):
        (tmp_path / source).write_bytes((RESULTS / source).read_bytes() + (b" " if source == name else b""))
    with pytest.raises(ValueError):
        module.load_context(tmp_path)

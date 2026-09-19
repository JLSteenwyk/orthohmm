from copy import deepcopy
from pathlib import Path

import pytest

from benchmark_tools import verify_frontier_overhead_provenance as module

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def fixture(index, panel="frontier_21838"):
    context = module.load_context(RESULTS, panel)
    spec = module.panel_spec(panel)
    plan, recipe = context["plan"], context["recipe"]
    task = plan["runs"][index]
    run = deepcopy(task["run"])
    run["gnu_time"] = dict(executable="/usr/bin/time", output=str(Path(run["measurement_directory"]).parent / "native.time.tsv"))
    order = plan["order"]
    originals = deepcopy(order["inputs_in_native_order"])
    copies = []
    if task["method"] == "orthofinder_full":
        for item in sorted(originals, key=lambda item: Path(item["path"]).name):
            copies.append(dict(item, path=str(Path(run["configuration"]["copy_inputs_to"]) / Path(item["path"]).name)))
    prepared = dict(run=run, measured_argv=module.time_command(run["native_argv"], run["gnu_time"]["output"]),
                    original_inputs=originals, copied_inputs=copies,
                    expected_native_basename_order=sorted(order["native_order"]) if copies else order["native_order"],
                    status="fresh_native_inputs_prepared", inference_started=False, preparation_wall_s=.1)
    runtime = [dict(item, records=count, scientific_execution_authorized=False, status="runtime_tree_identity_matches")
               for item, count in zip(plan["runtime_manifests"], (26673, 10066))]
    runtime.append(dict(path=str(module.ROOT / spec["recipe_manifest"]), sha256=spec["recipe_sha"],
                        records=len(recipe["records"]), scientific_execution_authorized=False, status="runtime_tree_identity_matches"))
    before = dict(runtime=runtime, original_inputs=originals, native_order=order["native_order"])
    after = deepcopy(before)
    if copies:
        after["copied_inputs"] = deepcopy(copies)
    collector = "measure_frontier_boundary_step.py" if task["mode"] == "boundary" else "measure_native_frontier_step.py"
    if spec.get("dual") and task["mode"] == "periodic":
        collector = "measure_native_dual_bracket_step.py"
    measured = dict(job_id=123, launched=["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
        str(module.ROOT / "envs/orthohmm/bin/python"), "-B",
        str(module.ROOT / spec["recipe_root"] / "benchmark_tools" / collector), "--worker", run["measurement_directory"]])
    wrapper = next(row["sha256"] for row in recipe["records"] if row["path"].endswith("/run_verified_slurm_measurement.py"))
    verified = dict(status="command_exited_zero", before=deepcopy(before), after=after,
                    measurement=deepcopy(measured), source_sha256=wrapper, scientific_results_admitted=False,
                    before_check_wall_s=1., after_check_wall_s=1.)
    receipt = dict(task=deepcopy(task), authorization=deepcopy(context["authorization"]),
                   authorization_sha256=spec["auth_sha"], plan_sha256=spec["plan_sha"], scientific_timings_admitted=False)
    if spec.get("dual"):
        receipt["recipe_sha256"] = spec["recipe_sha"]
    scheduler = (f"JobId=123 ArrayJobId={spec['array_id']} ArrayTaskId={index} JobState=COMPLETED ExitCode=0:0 Restarts=0 Requeue=0 "
                 "NodeList=spark-7ff0 OverSubscribe=NO MinMemoryNode=96G NumNodes=1 NumCPUs=20 CPUs/Task=20 "
                 f"Command={spec['scheduler_command']}")
    return [context, index, prepared, verified, receipt, measured, scheduler]


@pytest.mark.parametrize("index", range(18))
@pytest.mark.parametrize("panel", sorted(module.PANELS))
def test_all_tasks_bind_to_pinned_context(index, panel):
    args = fixture(index, panel)
    result = module.verify(*args)
    assert result["job_id"] == 123
    assert result["run"] == args[2]["run"]
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["receipt", "receipt_bool", "plan", "command", "input", "order", "runtime",
    "after", "wrapper", "measurement", "admitted", "duration", "worker", "scheduler_job", "scheduler_cpu",
    "scheduler_restart", "scheduler_exclusive", "scheduler_duplicate"])
@pytest.mark.parametrize("panel", sorted(module.PANELS))
def test_changed_provenance_rejected(fault, panel):
    args = fixture(0, panel)
    _, _, prepared, verified, receipt, measured, scheduler = args
    if fault == "receipt":
        receipt["authorization_sha256"] = "0"*64
    elif fault == "receipt_bool":
        receipt["authorization"]["execution_authorized"] = 1
    elif fault == "plan":
        receipt["task"]["mode"] = "periodic"
    elif fault == "command":
        prepared["run"]["native_argv"].append("--unexpected")
    elif fault == "input":
        prepared["original_inputs"][0]["sha256"] = "0"*64
    elif fault == "order":
        prepared["expected_native_basename_order"] = []
    elif fault == "runtime":
        verified["before"]["runtime"][0]["sha256"] = "0"*64
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
        args[-1] += " JobId=123"
    else:
        before, after = {"scheduler_job": ("JobId=123", "JobId=124"),
            "scheduler_cpu": ("CPUs/Task=20", "CPUs/Task=19"),
            "scheduler_restart": ("Restarts=0", "Restarts=1"),
            "scheduler_exclusive": ("OverSubscribe=NO", "OverSubscribe=YES")}[fault]
        args[-1] = scheduler.replace(before, after)
    with pytest.raises(ValueError):
        module.verify(*args)


def test_orthofinder_postcopy_proof_required():
    args = fixture(4)
    args[3]["after"].pop("copied_inputs")
    with pytest.raises(ValueError):
        module.verify(*args)


def test_retained_completed_scheduler_record_matches_submission():
    path = RESULTS / "dgx_frontier_overhead_scheduler_0_21838.txt"
    assert module.scheduler_identity(path.read_text(), 0) == 21839


@pytest.mark.parametrize("part", ["plan", "recipe", "authorization", "panel", "scheduler"])
def test_cross_panel_records_rejected(part):
    args = fixture(0, "pressure_21889")
    old = fixture(0)
    if part == "scheduler":
        args[-1] = old[-1]
    else:
        args[0][part] = old[0][part]
    with pytest.raises((ValueError, StopIteration)):
        module.verify(*args)


@pytest.mark.parametrize("key", ["plan_file", "recipe_file", "auth_file"])
def test_pressure_context_pins_input_bytes(tmp_path, key):
    spec = module.panel_spec("pressure_21889")
    for name in ("plan_file", "recipe_file", "auth_file"):
        data = (RESULTS / spec[name]).read_bytes()
        (tmp_path / spec[name]).write_bytes(data + (b" " if name == key else b""))
    with pytest.raises(ValueError):
        module.load_context(tmp_path, "pressure_21889")


def test_unregistered_panel_rejected():
    with pytest.raises(ValueError, match="Unknown"):
        module.load_context(RESULTS, "pressure_21869")

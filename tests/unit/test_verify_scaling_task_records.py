from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import verify_scaling_task_records as module

ROOT = Path(__file__).resolve().parents[2]
PLAN = ROOT / "benchmark_tools/results/dgx_root_context_scaling_plan_v2_20260920.json"


def save(path, value): path.write_text(json.dumps(value))


def fixture(tmp_path, index=0, code=0):
    plan, task, order = module.load_task(PLAN, index)
    directory = tmp_path / "run"
    (directory / "measurement").mkdir(parents=True)
    names = ("run_verified_slurm_measurement.py", "measure_scaling_root_context.py",
             "measure_root_context_scaling.py", "replay_scaling_root_context.py")
    entries = [dict(module.record(ROOT / "benchmark_tools" / name), kind="file",
                    path=str(module.RECIPE_ROOT / "benchmark_tools" / name)) for name in names]
    recipe = tmp_path / "recipe.json"
    save(recipe, dict(roots=[str(module.RECIPE_ROOT)], records=entries))
    sha = module.record(recipe)["sha256"]
    run = deepcopy(task["run"])
    remote = Path(run["measurement_directory"]).parent
    run["gnu_time"] = dict(executable="/usr/bin/time", output=str(remote / "native.time.tsv"))
    original = order["inputs_in_native_order"]
    copies = []
    if run["native_method"] == "orthofinder_full":
        copies = [dict(row, path=str(Path(run["configuration"]["copy_inputs_to"]) / Path(row["path"]).name))
                  for row in sorted(original, key=lambda r:Path(r["path"]).name)]
    prepared = dict(run=run, measured_argv=module.time_command(run["native_argv"], run["gnu_time"]["output"]),
        original_inputs=original, copied_inputs=copies, preparation_wall_s=1.,
        expected_native_basename_order=sorted(order["native_order"]) if copies else order["native_order"],
        status="fresh_native_inputs_prepared", inference_started=False)
    save(directory / "preparation.json", prepared)
    measured = dict(job_id=123, points=[{"root_context": {}}], native_wall_s=2.,
        native=dict(exit_code=code, timed_out=False, started_ns=1_000_000_000, finished_ns=3_000_000_000),
        status="command_exited_zero" if code == 0 else "command_failed",
        launched=["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
            plan["launcher_python"], "-B", str(module.RECIPE_ROOT / "benchmark_tools/measure_scaling_root_context.py"),
            "--worker", run["measurement_directory"]])
    runtime = [dict(row, records=count, scientific_execution_authorized=False, status="runtime_tree_identity_matches")
               for row,count in zip(plan["runtime_manifests"], (26673,10066))]
    runtime.append(dict(path=str(module.RECIPE_PATH), sha256=sha, records=len(entries),
                        scientific_execution_authorized=False, status="runtime_tree_identity_matches"))
    before = dict(runtime=runtime, original_inputs=original, native_order=order["native_order"])
    after = deepcopy(before)
    if copies: after["copied_inputs"] = copies
    verified = dict(status=measured["status"], before=before, after=after, measurement=measured,
        source_sha256=entries[0]["sha256"], scientific_results_admitted=False,
        before_check_wall_s=1., after_check_wall_s=1.)
    save(directory / "verification.json", verified)
    save(directory / "measurement/lineage_report.json", measured)
    save(directory / "measurement/root_context_report.json", dict(lineage_report=dict(
        module.record(directory / "measurement/lineage_report.json"), path="lineage_report.json")))
    return [PLAN, index, directory, recipe, sha, 123]


@pytest.mark.parametrize("index", range(27))
def test_all_original_tasks_bind_with_exact_orders(tmp_path, index):
    args = fixture(tmp_path, index)
    result = module.bind(*args)
    assert result["index"] == index
    assert result["run"]["dataset"]["proteomes"] == result["task"]["proteomes"]
    assert result["scientific_timings_admitted"] is False
    assert result["native_outputs_validated"] is False


@pytest.mark.parametrize("code", [7, -9, 124])
def test_failed_native_status_can_be_bound_without_becoming_success(tmp_path, code):
    result = module.bind(*fixture(tmp_path, 4, code))
    assert result["native_outcome"] == "exited_nonzero"
    assert result["measurement"]["native"]["exit_code"] == code


@pytest.mark.parametrize("fault", ["prepared_command", "input", "copies", "order", "runtime_count",
    "runtime_hash", "runtime_after", "wrapper", "worker", "embedded", "duration", "admission", "job",
    "context_hash", "recipe_root", "recipe_duplicate", "missing_source", "plan_hash"])
def test_contradictory_records_rejected(tmp_path, fault):
    args = fixture(tmp_path, 4)
    directory = args[2]
    if fault == "plan_hash": args[0] = args[3]
    elif fault.startswith("recipe_") or fault == "missing_source":
        data = json.loads(args[3].read_text())
        if fault == "recipe_root": data["roots"] = ["/other"]
        elif fault == "recipe_duplicate": data["records"].append(data["records"][0])
        else: data["records"].pop()
        save(args[3], data)
        args[4] = module.record(args[3])["sha256"]
    elif fault in {"prepared_command", "input", "copies", "order"}:
        path = directory / "preparation.json"
        data = json.loads(path.read_text())
        if fault == "prepared_command": data["measured_argv"].append("changed")
        elif fault == "input": data["original_inputs"][0]["sha256"] = "changed"
        elif fault == "copies": data["copied_inputs"] = []
        else: data["expected_native_basename_order"].reverse()
        save(path, data)
    elif fault == "context_hash":
        path = directory / "measurement/root_context_report.json"
        data = json.loads(path.read_text())
        data["lineage_report"]["sha256"] = "changed"
        save(path, data)
    else:
        path = directory / "verification.json"
        data = json.loads(path.read_text())
        if fault == "runtime_count": data["before"]["runtime"][0]["records"] += 1
        elif fault == "runtime_hash": data["before"]["runtime"][0]["sha256"] = "changed"
        elif fault == "runtime_after": data["after"]["copied_inputs"] = []
        elif fault == "wrapper": data["source_sha256"] = "changed"
        elif fault == "duration": data["before_check_wall_s"] = True
        elif fault == "admission": data["scientific_results_admitted"] = True
        elif fault == "embedded": data["measurement"]["native_wall_s"] = 3.
        elif fault in {"worker", "job"}:
            if fault == "worker": data["measurement"]["launched"][8] = "/wrong.py"
            else: data["measurement"]["job_id"] = True
            raw = directory / "measurement/lineage_report.json"
            save(raw, data["measurement"])
            save(directory / "measurement/root_context_report.json", dict(lineage_report=dict(module.record(raw), path="lineage_report.json")))
        save(path, data)
    with pytest.raises(ValueError): module.bind(*args)


def test_task_evidence_changed_during_binding_rejected(tmp_path, monkeypatch):
    args = fixture(tmp_path)
    original = module.native_outcome
    def changed(*values):
        result = original(*values)
        save(args[2] / "verification.json", {})
        return result
    monkeypatch.setattr(module, "native_outcome", changed)
    with pytest.raises(ValueError): module.bind(*args)


@pytest.mark.parametrize("name", ["preparation.json", "verification.json", "measurement/lineage_report.json"])
def test_direct_evidence_links_rejected(tmp_path, name):
    args = fixture(tmp_path)
    path = args[2] / name
    target = path.with_suffix(".retained")
    path.rename(target)
    path.symlink_to(target)
    with pytest.raises(ValueError, match="links"): module.bind(*args)

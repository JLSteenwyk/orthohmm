from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import verify_root_context_overhead_provenance as module
from benchmark_tools.prepare_frontier_overhead_panel import relocate
from tests.unit.test_verify_root_context_native_provenance import fixture as native_fixture, RESULTS


def fixture(tmp_path, index=0):
    plan = json.loads((RESULTS / "dgx_root_context_overhead_plan_20260919.json").read_text())
    task = plan["runs"][index]
    old_index = task["native_parent_index"]
    parent = tmp_path / "parent"
    parent.mkdir()
    old = native_fixture(parent, old_index)

    def transform(value):
        value = relocate(value, str(module.ROOT / "root_context_native_v1" / f"run_{old_index:02d}"),
                         str(module.ROOT / "root_context_overhead_v1" / f"run_{index:02d}"))
        return relocate(value, str(module.ROOT / "root_context_native_recipe_v1"), str(module.RECIPE_ROOT))

    recipe_path = tmp_path / "recipe.json"
    recipe_path.write_text(json.dumps(transform(old[0]["recipe"])))
    sha = hashlib.sha256(recipe_path.read_bytes()).hexdigest()
    context = module.load_context(RESULTS, recipe_path, sha)
    prep = transform(old[2])
    verified = transform(json.loads(old[3].read_text()))
    measured = transform(old[5])
    measured["points"] = [{"root_context": {}}] if task["arm"] == "root_context" else [{}]
    for phase in ("before", "after"):
        verified[phase]["runtime"][-1]["sha256"] = sha
        verified[phase]["runtime"][-1]["path"] = str(module.RECIPE_PATH)
    directory = tmp_path / "task"
    (directory / "measurement").mkdir(parents=True)
    (directory / "preparation.json").write_text(json.dumps(prep))
    verified["measurement"] = deepcopy(measured)
    (directory / "verification.json").write_text(json.dumps(verified))
    (directory / "measurement/lineage_report.json").write_text(json.dumps(measured))
    job = old[-1]
    scheduler = old[-2].replace(str(module.ROOT / "root_context_native_recipe_v1"), str(module.RECIPE_ROOT)).replace("TimeLimit=01:00:00", "TimeLimit=05:00:00").replace(
        "run_dgx_root_context_native.sh", "run_dgx_root_context_overhead.sh") + " Partition=spark"
    receipt = {key: task[key] for key in ("index", "block", "pair", "arm", "method")}
    receipt.update(task=deepcopy(task), job_id=job, scientific_timings_admitted=False,
        status="measurement_completed", wrapper_status="command_exited_zero", native_wall_s=15.)
    args = [context, index, directory, receipt, scheduler, job]
    refresh(args)
    return args


def refresh(args):
    context, index, directory, receipt, *_ = args
    task = context["plan"]["runs"][index]
    remote = Path(task["run"]["measurement_directory"])
    for key, path, target in (("verification", directory / "verification.json", remote.parent / "verification.json"),
                              ("lineage_report", directory / "measurement/lineage_report.json", remote / "lineage_report.json")):
        receipt[key] = dict(module.record(path), path=str(target))
    if task["arm"] == "root_context":
        root = directory / "measurement/root_context_report.json"
        root.write_text(json.dumps(dict(lineage_report=dict(receipt["lineage_report"], path="lineage_report.json"))))
        receipt["root_context_report"] = dict(module.record(root), path=str(remote / "root_context_report.json"))


@pytest.mark.parametrize("index", range(18))
def test_all_prescribed_task_and_arm_identities_bind(tmp_path, index):
    args = fixture(tmp_path, index)
    result = module.verify(*args)
    assert result["index"] == index and result["arm"] == args[3]["arm"]
    assert result["scientific_timings_admitted"] is False
    assert len(result["evidence"]) == (4 if result["arm"] == "root_context" else 3)


@pytest.mark.parametrize("fault", ["pair", "arm", "task", "receipt_bool", "receipt_hash", "extra_root",
    "lineage_bytes", "root_bytes", "root_link", "point_arm", "prepared", "input", "copies",
    "runtime", "embedded", "wrapper", "worker", "exit", "duration", "scheduler_limit",
    "scheduler_duplicate", "scheduler_partition", "scheduler_running", "scheduler_array"])
def test_tampered_provenance_rejected(tmp_path, fault):
    args = fixture(tmp_path, 4 if fault == "extra_root" else 5)
    _, _, directory, receipt, _, _ = args
    if fault in {"pair", "arm", "task", "receipt_bool", "receipt_hash"}:
        if fault == "pair":
            receipt["pair"] += 1
        elif fault == "arm":
            receipt["arm"] = "lineage"
        elif fault == "task":
            receipt["task"]["index"] = 0
        elif fault == "receipt_bool":
            receipt["scientific_timings_admitted"] = 0
        else:
            receipt["lineage_report"]["sha256"] = "0"*64
    elif fault == "extra_root":
        (directory / "measurement/root_context_report.json").write_text("{}")
    elif fault in {"lineage_bytes", "root_bytes"}:
        path = directory / "measurement" / ("lineage_report.json" if fault == "lineage_bytes" else "root_context_report.json")
        path.write_bytes(path.read_bytes()+b" ")
    elif fault == "root_link":
        path = directory / "measurement/root_context_report.json"
        value = json.loads(path.read_text())
        value["lineage_report"]["sha256"] = "0"*64
        path.write_text(json.dumps(value))
        receipt["root_context_report"].update(bytes=path.stat().st_size, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    elif fault in {"prepared", "input", "copies"}:
        path = directory / "preparation.json"
        value = json.loads(path.read_text())
        if fault == "prepared":
            value["measured_argv"].append("wrong")
        elif fault == "input":
            value["original_inputs"][0]["sha256"] = "0"*64
        else:
            value["copied_inputs"] = []
        path.write_text(json.dumps(value))
    elif fault in {"runtime", "embedded", "wrapper", "duration", "worker", "exit", "point_arm"}:
        path = directory / "verification.json"
        value = json.loads(path.read_text())
        if fault == "runtime":
            value["after"]["runtime"][0]["records"] += 1
        elif fault == "embedded":
            value["measurement"]["job_id"] += 1
        elif fault == "wrapper":
            value["source_sha256"] = "0"*64
        elif fault == "duration":
            value["before_check_wall_s"] = -1
        else:
            if fault == "worker":
                value["measurement"]["launched"][8] = "/wrong.py"
            elif fault == "exit":
                value["measurement"]["native"]["exit_code"] = False
            else:
                value["measurement"]["points"] = [{}]
            (directory / "measurement/lineage_report.json").write_text(json.dumps(value["measurement"]))
        path.write_text(json.dumps(value))
        refresh(args)
    elif fault == "scheduler_duplicate":
        args[-2] += " JobId=12345"
    elif fault == "scheduler_array":
        args[-2] += " ArrayJobId=12345 ArrayTaskId=0"
    else:
        old, new = {"scheduler_limit": ("TimeLimit=05:00:00", "TimeLimit=01:00:00"),
            "scheduler_partition": ("Partition=spark", "Partition=other"),
            "scheduler_running": ("JobState=COMPLETED", "JobState=RUNNING")}[fault]
        args[-2] = args[-2].replace(old, new)
    with pytest.raises(ValueError):
        module.verify(*args)


@pytest.mark.parametrize("index", [-1, 18, True, 1.0])
def test_invalid_task_indices(tmp_path, index):
    args = fixture(tmp_path)
    args[1] = index
    with pytest.raises(ValueError):
        module.verify(*args)


def test_completed_task_in_later_failed_allocation(tmp_path):
    args = fixture(tmp_path)
    args[-2] = args[-2].replace("JobState=COMPLETED", "JobState=FAILED").replace("ExitCode=0:0", "ExitCode=1:0")
    with pytest.raises(ValueError):
        module.verify(*args)
    assert module.verify(*args, require_completed=False)["index"] == 0


@pytest.mark.parametrize("changed", ["dgx_root_context_overhead_plan_20260919.json",
    "dgx_root_context_native_plan_20260919.json", "ROOT_CONTEXT_OVERHEAD_PROTOCOL_20260919.md"])
def test_frozen_derivation_rejects_source_changes(tmp_path, changed):
    args = fixture(tmp_path)
    for name in ("dgx_root_context_overhead_plan_20260919.json", "dgx_root_context_native_plan_20260919.json",
                 "ROOT_CONTEXT_OVERHEAD_PROTOCOL_20260919.md"):
        (tmp_path / name).write_bytes((RESULTS / name).read_bytes() + (b" " if name == changed else b""))
    with pytest.raises(ValueError):
        module.load_context(tmp_path, tmp_path / "recipe.json", args[0]["recipe_sha256"])


@pytest.mark.parametrize("fault", ["bytes", "duplicate"])
def test_recipe_integrity(tmp_path, fault):
    args = fixture(tmp_path)
    path = tmp_path / "recipe.json"
    sha = args[0]["recipe_sha256"]
    if fault == "bytes":
        path.write_bytes(path.read_bytes()+b" ")
    else:
        recipe = json.loads(path.read_text())
        recipe["records"].append(recipe["records"][0])
        path.write_text(json.dumps(recipe))
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError):
        module.load_context(RESULTS, path, sha)


@pytest.mark.parametrize("fault", ["missing", "symlink"])
def test_missing_or_symlinked_supplementary_report(tmp_path, fault):
    args = fixture(tmp_path, 1)
    path = args[2] / "measurement/root_context_report.json"
    target = tmp_path / "report.json"
    path.rename(target)
    if fault == "symlink":
        path.symlink_to(target)
    with pytest.raises((ValueError, OSError)):
        module.verify(*args)

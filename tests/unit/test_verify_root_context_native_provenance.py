from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import verify_root_context_native_provenance as module
from tests.unit.test_verify_lineage_native_provenance import fixture as lineage_fixture

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def fixture(tmp_path, index=0):
    old = lineage_fixture(index)

    def relocate(value):
        if isinstance(value, str):
            return value.replace("lineage_native_v1", "root_context_native_v1").replace(
                "lineage_native_recipe_v1", "root_context_native_recipe_v1")
        if isinstance(value, dict):
            return {k: relocate(v) for k, v in value.items()}
        if isinstance(value, list):
            return [relocate(v) for v in value]
        return value

    recipe = relocate(old[0]["recipe"])
    recipe_path = tmp_path / "recipe.json"
    recipe_path.write_text(json.dumps(recipe))
    sha = hashlib.sha256(recipe_path.read_bytes()).hexdigest()
    context = module.load_context(RESULTS, recipe_path, sha)
    prepared, verified, _, measured = relocate(old[2:6])
    job = 12345
    measured.update(job_id=job, native_wall_s=15., native=dict(exit_code=0, timed_out=False))
    verified["measurement"] = deepcopy(measured)
    for phase in ("before", "after"):
        verified[phase]["runtime"][-1]["sha256"] = sha
    path = tmp_path / "verification.json"
    task = context["plan"]["runs"][index]
    receipt = dict(index=index, task=deepcopy(task), job_id=job, scientific_timings_admitted=False,
        status="measurement_completed", wrapper_status="command_exited_zero", native_wall_s=15.)
    scheduler = (f"JobId={job} JobState=COMPLETED ExitCode=0:0 Restarts=0 Requeue=0 "
        "NodeList=spark-7ff0 OverSubscribe=NO MinMemoryNode=96G NumNodes=1 NumCPUs=20 "
        f"NumTasks=1 CPUs/Task=20 TimeLimit=01:00:00 Command={module.SUBMISSION_SCRIPT} "
        f"WorkDir={module.RECIPE_ROOT}")
    args = [context, index, prepared, path, receipt, measured, scheduler, job]
    save_verification(args, verified)
    return args


def save_verification(args, verified):
    args[3].write_text(json.dumps(verified))
    raw = args[3].read_bytes()
    args[4]["verification"] = dict(
        path=str(Path(args[0]["plan"]["runs"][args[1]]["run"]["measurement_directory"]).parent / "verification.json"),
        bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


@pytest.mark.parametrize("index", range(3))
def test_native_tasks_bind_to_one_job(tmp_path, index):
    args = fixture(tmp_path, index)
    result = module.verify(*args)
    assert result["job_id"] == 12345
    assert result["run"] == args[2]["run"]
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["receipt_hash", "receipt_path", "receipt_bool", "receipt_task",
    "file_bytes", "command", "input", "order", "copies", "runtime", "wrapper", "worker",
    "embedded", "duration", "native_failure", "native_timeout", "native_bool", "job_bool",
    "unrun", "duplicate_scheduler", "array", "resources", "time_limit", "script", "cwd"])
def test_tampering_is_rejected(tmp_path, fault):
    args = fixture(tmp_path, 2)
    _, _, prep, path, receipt, measured, _, _ = args
    verified = json.loads(path.read_text())
    if fault == "receipt_hash":
        receipt["verification"]["sha256"] = "0" * 64
    elif fault == "receipt_path":
        receipt["verification"]["path"] = str(path)
    elif fault == "receipt_bool":
        receipt["scientific_timings_admitted"] = 0
    elif fault == "receipt_task":
        receipt["task"]["index"] = 0
    elif fault == "file_bytes":
        path.write_bytes(path.read_bytes() + b" ")
    elif fault == "command":
        prep["run"]["native_argv"].append("--wrong")
    elif fault == "input":
        prep["original_inputs"][0]["sha256"] = "0" * 64
    elif fault == "order":
        prep["expected_native_basename_order"] = []
    elif fault == "unrun":
        receipt["status"] = "not_run_after_failure"
    elif fault == "job_bool":
        args[-1] = True
    elif fault in ("copies", "runtime", "wrapper", "worker", "embedded", "duration",
                    "native_failure", "native_timeout", "native_bool"):
        if fault == "copies":
            verified["after"].pop("copied_inputs")
        elif fault == "runtime":
            verified["before"]["runtime"][0]["records"] += 1
        elif fault == "wrapper":
            verified["source_sha256"] = "0" * 64
        elif fault == "embedded":
            verified["measurement"]["job_id"] += 1
        elif fault == "duration":
            verified["before_check_wall_s"] = -1
        else:
            if fault == "worker":
                measured["launched"][8] = "/wrong/worker.py"
            elif fault == "native_failure":
                measured["native"]["exit_code"] = 1
            elif fault == "native_bool":
                measured["native"]["exit_code"] = False
            else:
                measured["native"]["timed_out"] = True
            verified["measurement"] = deepcopy(measured)
        save_verification(args, verified)
    elif fault == "duplicate_scheduler":
        args[-2] += " JobId=12345"
    elif fault == "array":
        args[-2] += " ArrayJobId=12345 ArrayTaskId=2"
    else:
        before, after = {"resources": ("NumCPUs=20", "NumCPUs=19"),
            "time_limit": ("TimeLimit=01:00:00", "TimeLimit=00:15:00"),
            "script": (f"Command={module.SUBMISSION_SCRIPT}", "Command=/wrong.sh"),
            "cwd": (f"WorkDir={module.RECIPE_ROOT}", "WorkDir=/tmp")}[fault]
        args[-2] = args[-2].replace(before, after)
    with pytest.raises(ValueError):
        module.verify(*args)


@pytest.mark.parametrize("index", [-1, 3, True, 1.0])
def test_invalid_task_index(tmp_path, index):
    args = fixture(tmp_path)
    args[1] = index
    with pytest.raises(ValueError):
        module.verify(*args)


def test_recipe_byte_drift(tmp_path):
    args = fixture(tmp_path)
    recipe = tmp_path / "recipe.json"
    recipe.write_bytes(recipe.read_bytes() + b" ")
    with pytest.raises(ValueError):
        module.load_context(RESULTS, recipe, args[0]["recipe_sha256"])


def test_duplicate_recipe_record(tmp_path):
    args = fixture(tmp_path)
    recipe = deepcopy(args[0]["recipe"])
    recipe["records"].append(deepcopy(recipe["records"][0]))
    path = tmp_path / "duplicate.json"
    path.write_text(json.dumps(recipe))
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="Duplicate"):
        module.load_context(RESULTS, path, sha)


@pytest.mark.parametrize("changed", ["dgx_root_context_native_plan_20260919.json",
    "dgx_lineage_native_plan_20260919.json", "ROOT_CONTEXT_NATIVE_PROTOCOL_20260919.md"])
def test_frozen_plan_inputs_reject_byte_changes(tmp_path, changed):
    args = fixture(tmp_path)
    for name in ("dgx_root_context_native_plan_20260919.json",
                 "dgx_lineage_native_plan_20260919.json", "ROOT_CONTEXT_NATIVE_PROTOCOL_20260919.md"):
        (tmp_path / name).write_bytes((RESULTS / name).read_bytes() + (b" " if name == changed else b""))
    with pytest.raises(ValueError):
        module.load_context(tmp_path, tmp_path / "recipe.json", args[0]["recipe_sha256"])

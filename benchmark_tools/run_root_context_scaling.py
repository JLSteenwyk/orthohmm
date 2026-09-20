"""Execute one explicitly authorized replacement-scaling task; never submit or retry."""

import argparse
import hashlib
from pathlib import Path
import subprocess
import sys
import time

from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.measure_root_context_scaling import load_task, measure_task, PLAN_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_root_context_scaling import OUTPUT_ROOT
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_root_context_native import preflight
from benchmark_tools.verify_scaling_allocation import validate, SUBMISSION_SCRIPT
from benchmark_tools.verify_scaling_task_records import RECIPE_ROOT, RECIPE_PATH


def authorization(path, sha, recipe_sha, index, job):
    """Bind external authorization and preflight receipts, not their scientific validity."""
    documents = [record(path)]
    permit = read_pinned(path, sha)
    expected = dict(schema="root_context_scaling_authorization_v1", execution_authorized=True,
        plan_sha256=PLAN_SHA, recipe_sha256=recipe_sha, index=index, job_id=job)
    if any(type(permit.get(k)) is not type(v) or permit[k] != v for k, v in expected.items()):
        raise ValueError("Explicit per-job execution authorization differs")
    policy_ref, ready_ref = permit["environment_policy"], permit["environment_preflight"]
    documents.extend(record(ref["path"]) for ref in (policy_ref, ready_ref))
    policy = read_pinned(policy_ref["path"], policy_ref["sha256"])
    ready = read_pinned(ready_ref["path"], ready_ref["sha256"])
    if (policy.get("status") != "approved_frozen_environment_policy"
            or not isinstance(policy.get("approval_reference"), str)
            or not policy["approval_reference"].strip()):
        raise ValueError("Require separately approved frozen environmental policy")
    expected_ready = dict(status="environment_preflight_passed", job_id=job,
        policy_sha256=policy_ref["sha256"], whole_run_observer_ready=True)
    if any(type(ready.get(k)) is not type(v) or ready[k] != v for k, v in expected_ready.items()):
        raise ValueError("Environmental preflight or whole-run observer receipt differs")
    for item in documents:
        check(item)
    return documents


def select(plan_path, index, recipe_path, recipe_sha):
    plan, task, _ = load_task(plan_path, index)
    if Path(__file__).resolve().parent.parent != RECIPE_ROOT or Path(recipe_path).resolve() != RECIPE_PATH:
        raise ValueError("Require frozen deployed recipe location")
    recipe = read_pinned(recipe_path, recipe_sha)
    if recipe["roots"] != [str(RECIPE_ROOT)]:
        raise ValueError("Recipe root differs")
    rows = recipe["records"]
    if len({row["path"] for row in rows}) != len(rows):
        raise ValueError("Duplicate recipe records")
    files = {row["path"]: row for row in rows if row["kind"] == "file"}
    required = [*Path(__file__).resolve().parent.glob("*.py"), Path(plan_path).resolve(), SUBMISSION_SCRIPT]
    for path in required:
        data = path.read_bytes()
        item = files.get(str(path), {})
        if item.get("sha256") != hashlib.sha256(data).hexdigest() or item.get("bytes") != len(data):
            raise ValueError("Executor source, script or plan not pinned in recipe")
    return plan, task


def execute(plan_path, index, recipe_path, recipe_sha, authorization_path, authorization_sha):
    plan_path, recipe_path, authorization_path = map(Path, (plan_path, recipe_path, authorization_path))
    plan, task = select(plan_path, index, recipe_path, recipe_sha)
    job = preflight(plan)
    if Path(sys.pycache_prefix).is_symlink():
        raise ValueError("Observer cache prefix must not be a dangling link")
    documents = [record(plan_path), record(recipe_path), *authorization(
        authorization_path, authorization_sha, recipe_sha, index, job)]
    receipts = OUTPUT_ROOT / "sessions" / f"task_{index:02d}"
    receipts.mkdir(parents=True, exist_ok=False)
    result = dict(status="prelaunch_pending", index=index, job_id=job, task=task,
        plan_sha256=PLAN_SHA, recipe_sha256=recipe_sha, evidence=documents,
        scientific_timings_admitted=False, native_outputs_validated=False,
        scheduler_terminal_verified=False, environmental_validity_established=False,
        next_submission_authorized=False, automatic_retry=False)
    try:
        command = ["scontrol", "show", "job", str(job), "--oneliner"]
        started = time.time_ns()
        controller = subprocess.run(command, capture_output=True, text=True, timeout=15)
        save(receipts / "controller_before.json", dict(command=command,
            started_unix_ns=started, finished_unix_ns=time.time_ns(),
            returncode=controller.returncode, stdout=controller.stdout, stderr=controller.stderr))
        if controller.returncode:
            raise ValueError("Controller preflight failed")
        result["allocation"] = validate(controller.stdout, job, "running")
        for item in documents:
            check(item)
        save(receipts / "launch.json", dict(index=index, job_id=job, started_unix_ns=time.time_ns(),
            evidence=documents, source=record(__file__), scientific_timings_admitted=False))
        measured = measure_task(plan_path, index, recipe_path, recipe_sha, job)
        result["wrapper_status"] = measured["status"]
        result["verification"] = record(Path(task["run"]["measurement_directory"]).parent / "verification.json")
        for item in documents:
            check(item)
        select(plan_path, index, recipe_path, recipe_sha)
        result["status"] = "measurement_returned_pending_audit"
    except BaseException as error:
        result.update(status="executor_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        result["finished_unix_ns"] = time.time_ns()
        save(receipts / "result.json", result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "recipe", "authorization"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--recipe-sha", required=True)
    parser.add_argument("--authorization-sha", required=True)
    args = parser.parse_args(argv)
    result = execute(args.plan, args.index, args.recipe, args.recipe_sha,
                     args.authorization, args.authorization_sha)
    return 0 if result["wrapper_status"] == "command_exited_zero" else 1


if __name__ == "__main__":
    raise SystemExit(main())

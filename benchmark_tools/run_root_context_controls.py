"""Run the frozen root-context panel with no retries or timing admission."""

import argparse
import hashlib
import os
from pathlib import Path
import platform
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.root_context_control_design import ORDER, check_service_scope, unit_name
from benchmark_tools.run_full_node_controls import preflight as runtime_preflight, PLAN_SHA
from benchmark_tools.run_root_context_control_trial import trial
from benchmark_tools.run_verified_slurm_measurement import run_checked

PROTOCOL_SHA = "4fe2ae15ddd11337cc9537b0042eb74693cd636bffb0ba2d75be5148594d5887"


def summarize(rows):
    expected = [(block, mode) for block, modes in enumerate(ORDER) for mode in modes]
    if [(row["block"], row["mode"]) for row in rows] != expected:
        raise ValueError("Root-control inventory/order differs")
    if [row["index"] for row in rows] != list(range(12)):
        raise ValueError("Root-control indices differ")
    valid = all(row["status"] == "root_context_workload_validated" for row in rows)
    responses = [dict(index=row["index"], detected=row.get("positive_control_response")
        if row["status"] == "root_context_workload_validated" else None)
        for row in rows if row["mode"] == "user-contended"]
    return dict(all_workloads_valid=valid, positive_controls=responses,
        all_positive_controls_detected=all(row["detected"] is True for row in responses),
        scientific_timings_admitted=False)


def panel(output, job, manager):
    check_service_scope("0::" + manager + "/" + unit_name(job, 0), manager, unit_name(job, 0))
    output.mkdir(exist_ok=False)
    rows, failed_index = [], None
    for block, modes in enumerate(ORDER):
        for mode in modes:
            index = len(rows)
            directory = output / f"trial_{index:02d}"
            identity = dict(block=block, mode=mode, index=index, job_id=job,
                            scientific_timings_admitted=False)
            if failed_index is not None:
                row = dict(identity, status="not_run_after_failure", failed_index=failed_index)
            else:
                try:
                    result = trial(directory, mode, job, index, manager)
                    if result.get("status") != "root_context_workload_validated":
                        raise ValueError("Unexpected trial status")
                    row = dict(result, **identity)
                except Exception as error:
                    failed_index = index
                    row = dict(identity, status="failed", error_type=type(error).__name__, error=str(error))
            directory.mkdir(exist_ok=True)
            save(directory / "panel_trial.json", row)
            rows.append(row)
            save(output / f"progress_{index:02d}.json", dict(trials=rows, job_id=job, failed_index=failed_index))
    result = dict(status="root_context_controls_completed" if failed_index is None else "root_context_controls_stopped",
        trials=rows, job_id=job, manager=manager, summary=summarize(rows),
        protocol_sha256=PROTOCOL_SHA, scientific_timings_admitted=False, publication_ready=False)
    save(output / "result.json", result)
    return result


def preflight(recipe_path, recipe_sha):
    # Reuse environment checks only; no native-method commands from the old plan run here.
    plan, job = runtime_preflight(recipe_path, recipe_sha)
    recipe = read_pinned(recipe_path, recipe_sha)
    files = {row["path"]: row for row in recipe["records"] if row["kind"] == "file"}
    base = Path(__file__).resolve().parent
    protocol = base / "results/ROOT_CPU_CONTEXT_CONTROL_PROTOCOL_20260919.md"
    for path in (protocol, base / "run_dgx_root_context_controls.sh"):
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        if files.get(str(path), {}).get("sha256") != sha:
            raise ValueError("Root protocol or launch script not pinned")
        if path == protocol and sha != PROTOCOL_SHA:
            raise ValueError("Frozen root protocol changed")
    if os.environ.get("PYTHONNOUSERSITE") != "1":
        raise ValueError("Require disabled user site packages")
    manager = subprocess.check_output(["systemctl", "--user", "show", "--property=ControlGroup", "--value"],
                                     text=True, timeout=10).strip()
    check_service_scope("0::" + manager + "/" + unit_name(job, 0), manager, unit_name(job, 0))
    if not (Path("/sys/fs/cgroup") / manager.lstrip("/")).is_dir():
        raise ValueError("User-manager scope unavailable")
    return plan, job, manager


def run(output, recipe_path, recipe_sha):
    plan, job, manager = preflight(recipe_path, recipe_sha)
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    os.environ["PYTHONPYCACHEPREFIX"] = sys.pycache_prefix
    specifications = [(row["path"], row["sha256"]) for row in plan["runtime_manifests"]]
    specifications.append((str(recipe_path), recipe_sha))
    result = run_checked(specifications, output, lambda directory: panel(directory, job, manager))
    save(output / "launch.json", dict(job_id=job, manager=manager, recipe_sha256=recipe_sha,
        protocol_sha256=PROTOCOL_SHA, runtime_plan_sha256=PLAN_SHA, executable=sys.executable,
        python=sys.version, kernel=platform.release(), host=platform.node(),
        source_directory=str(Path(__file__).resolve().parent.parent),
        scientific_timings_admitted=False))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--recipe-sha", required=True)
    args = parser.parse_args()
    result = run(args.output.resolve(), args.recipe.resolve(), args.recipe_sha)
    summary = result.get("measurement", {}).get("summary", {})
    raise SystemExit(0 if result["status"] == "root_context_controls_completed"
        and summary.get("all_workloads_valid") is True
        and summary.get("all_positive_controls_detected") is True else 1)

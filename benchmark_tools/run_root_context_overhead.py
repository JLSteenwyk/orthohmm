"""Execute the frozen 18-task incremental collector panel serially."""

import argparse
import hashlib
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_tools.launch_dgx_native_run import read_pinned, native_enumerator
from benchmark_tools.measure_native_lineage_step import measure as measure_lineage
from benchmark_tools.measure_native_root_context import measure_native_run as measure_root
from benchmark_tools.measure_native_scaling_run import measure_run
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_root_context_overhead import build, OUTPUT_ROOT, PROTOCOL_SHA
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_root_context_native import preflight
from benchmark_tools.verify_lineage_native_provenance import same

PLAN_SHA = "02475e4a4290664f776d430b041bd65ccfa640221d07f95b1709b53678a91873"
COLLECTORS = {"lineage": measure_lineage, "root_context": measure_root}


def select(plan_path, recipe_path, recipe_sha):
    folder = Path(__file__).resolve().parent
    plan = read_pinned(plan_path, PLAN_SHA)
    expected = build(folder / "results/dgx_root_context_native_plan_20260919.json",
                     folder / "results/ROOT_CONTEXT_OVERHEAD_PROTOCOL_20260919.md")
    if not same(plan, expected):
        raise ValueError("Overhead plan differs from frozen derivation")
    recipe = read_pinned(recipe_path, recipe_sha)
    paths = [row["path"] for row in recipe["records"]]
    if len(paths) != len(set(paths)):
        raise ValueError("Duplicate recipe record")
    files = {row["path"]: row for row in recipe["records"] if row["kind"] == "file"}
    required = [*folder.glob("*.py"), plan_path.resolve(), folder / "run_dgx_root_context_overhead.sh",
        folder / "results/dgx_root_context_native_plan_20260919.json",
        folder / "results/ROOT_CONTEXT_OVERHEAD_PROTOCOL_20260919.md"]
    for path in required:
        if files.get(str(path), {}).get("sha256") != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError("Overhead source or input not pinned")
    return plan


def execute_task(plan, task, recipe_path, recipe_sha, job):
    collector = COLLECTORS[task["arm"]]
    run = task["run"]
    cache = OUTPUT_ROOT / f"cache_{task['index']:02d}"
    if cache.exists() or cache.is_symlink():
        raise ValueError("Task cache prefix must be absent")
    overrides = dict(plan["environment_overrides"],
        PATH=os.pathsep.join(plan["environment_paths"][run["environment_role"]]),
        PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(cache))
    keys = set(overrides) | set(plan["unset_environment"])
    previous, cwd = {key: os.environ.get(key) for key in keys}, Path.cwd()
    try:
        for key in plan["unset_environment"]:
            os.environ.pop(key, None)
        os.environ.update(overrides)
        os.chdir(run["cwd"])
        specifications = [(row["path"], row["sha256"]) for row in plan["runtime_manifests"]]
        specifications.append((str(recipe_path), recipe_sha))
        return measure_run(run, plan["order"], specifications, native_enumerator(plan["enumerator"]),
                           collector, job, timeout_s=900)
    finally:
        os.chdir(cwd)
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def panel(plan_path, recipe_path, recipe_sha):
    plan = select(plan_path, recipe_path, recipe_sha)
    job = preflight(plan)
    output = OUTPUT_ROOT
    output.mkdir(exist_ok=False)
    save(output / "launch.json", dict(job_id=job, plan_sha256=PLAN_SHA, recipe_sha256=recipe_sha,
        protocol_sha256=PROTOCOL_SHA, executable=sys.executable, python=sys.version,
        source_directory=str(Path(__file__).resolve().parent.parent), uname=list(os.uname()),
        affinity=sorted(os.sched_getaffinity(0)), scientific_timings_admitted=False))
    rows, failed = [], None
    for task in plan["runs"]:
        index = task["index"]
        row = {key: task[key] for key in ("index", "block", "pair", "arm", "method")}
        row.update(task=task, job_id=job, scientific_timings_admitted=False)
        if failed is not None:
            row.update(status="not_run_after_failure", failed_index=failed)
        else:
            try:
                if not same(select(plan_path, recipe_path, recipe_sha), plan):
                    raise ValueError("Plan changed before task")
                result = execute_task(plan, task, recipe_path, recipe_sha, job)
                row["wrapper_status"] = result["status"]
                directory = Path(task["run"]["measurement_directory"])
                row["verification"] = record(directory.parent / "verification.json")
                if result["status"] != "command_exited_zero":
                    raise ValueError("Native preparation, measurement or runtime verification failed")
                native = result["measurement"]["native"]
                if type(native["exit_code"]) is not int or native["exit_code"] != 0 or native["timed_out"] is not False:
                    raise ValueError("Native command failed or timed out")
                row["lineage_report"] = record(directory / "lineage_report.json")
                supplementary = directory / "root_context_report.json"
                if task["arm"] == "root_context":
                    row["root_context_report"] = record(supplementary)
                elif supplementary.exists():
                    raise ValueError("Unexpected supplementary report in lineage arm")
                if not same(select(plan_path, recipe_path, recipe_sha), plan):
                    raise ValueError("Plan changed after task")
                row.update(status="measurement_completed", native_wall_s=result["measurement"]["native_wall_s"])
            except Exception as error:
                failed = index
                row.update(status="failed", error_type=type(error).__name__, error=str(error))
        save(output / f"task_{index:02d}.json", row)
        rows.append(row)
        save(output / f"progress_{index:02d}.json", dict(tasks=rows, failed_index=failed, job_id=job))
    result = dict(status="root_context_overhead_completed" if failed is None else "root_context_overhead_stopped",
        tasks=rows, job_id=job, plan_sha256=PLAN_SHA, recipe_sha256=recipe_sha,
        scientific_timings_admitted=False, native_outputs_validated=False, publication_ready=False)
    save(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--recipe-sha", required=True)
    args = parser.parse_args()
    result = panel(args.plan.resolve(), args.recipe.resolve(), args.recipe_sha)
    raise SystemExit(0 if result["status"] == "root_context_overhead_completed" else 1)

"""Execute the three frozen native integration steps serially, retaining failures."""

import argparse
import hashlib
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_tools.launch_dgx_native_run import read_pinned, native_enumerator
from benchmark_tools.measure_native_root_context import measure_native_run
from benchmark_tools.measure_native_scaling_run import measure_run
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_root_context_native import ROOT, build, PROTOCOL_SHA
from benchmark_tools.probe_dgx_step_separation import save

PLAN_SHA = "9b99cd810aaf7e040cda0242dd9b6d1da82bcb231c22657b8dc0108c4c4178ad"


def select(plan_path, recipe_path, recipe_sha):
    folder = Path(__file__).resolve().parent
    plan = read_pinned(plan_path, PLAN_SHA)
    expected = build(folder / "results/dgx_lineage_native_plan_20260919.json",
                     folder / "results/ROOT_CONTEXT_NATIVE_PROTOCOL_20260919.md")
    if plan != expected:
        raise ValueError("Native root-context plan differs from frozen derivation")
    recipe = read_pinned(recipe_path, recipe_sha)
    files = {r["path"]: r for r in recipe["records"] if r["kind"] == "file"}
    if len(files) != sum(r["kind"] == "file" for r in recipe["records"]):
        raise ValueError("Duplicate recipe file")
    required = [*folder.glob("*.py"), plan_path.resolve(), folder / "run_dgx_root_context_native.sh",
                folder / "results/dgx_lineage_native_plan_20260919.json",
                folder / "results/ROOT_CONTEXT_NATIVE_PROTOCOL_20260919.md"]
    for path in required:
        if files.get(str(path), {}).get("sha256") != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError("Native launcher source or input not pinned")
    return plan


def preflight(plan):
    if (Path(sys.executable).absolute() != Path(plan["launcher_python"])
            or os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "20"
            or os.environ.get("SLURM_MEM_PER_NODE") != "98304" or len(os.sched_getaffinity(0)) != 20):
        raise ValueError("Require pinned DGX20CPU/96GiB interpreter and allocation")
    if (not sys.dont_write_bytecode or os.environ.get("PYTHONHASHSEED") != "0"
            or os.environ.get("PYTHONNOUSERSITE") != "1" or not sys.pycache_prefix
            or Path(sys.pycache_prefix).exists()):
        raise ValueError("Require isolated absent cache, disabled writes/user-site and fixed hash seed")
    if any(os.environ.get(k) for k in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT", "PYTHONPATH")) or Path("/etc/ld.so.preload").exists():
        raise ValueError("Unexpected loader or Python override")
    job = int(os.environ["SLURM_JOB_ID"])
    if job <= 0:
        raise ValueError("Invalid job identity")
    return job


def execute_task(plan, task, recipe_path, recipe_sha, job):
    run = task["run"]
    cache = ROOT / "root_context_native_v1" / f"cache_{task['index']}"
    if cache.exists():
        raise ValueError("Native cache prefix must be absent")
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
        specifications = [(r["path"], r["sha256"]) for r in plan["runtime_manifests"]]
        specifications.append((str(recipe_path), recipe_sha))
        return measure_run(run, plan["order"], specifications, native_enumerator(plan["enumerator"]),
                           measure_native_run, job, timeout_s=900)
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
    output = ROOT / "root_context_native_v1"
    output.mkdir(exist_ok=False)
    save(output / "launch.json", dict(job_id=job, plan_sha256=PLAN_SHA, recipe_sha256=recipe_sha,
        protocol_sha256=PROTOCOL_SHA, executable=sys.executable, python=sys.version,
        source_directory=str(Path(__file__).resolve().parent.parent),
        uname=list(os.uname()), affinity=sorted(os.sched_getaffinity(0)), scientific_timings_admitted=False))
    rows, failed = [], None
    for task in plan["runs"]:
        index = task["index"]
        row = dict(index=index, task=task, job_id=job, scientific_timings_admitted=False)
        if failed is not None:
            row.update(status="not_run_after_failure", failed_index=failed)
        else:
            try:
                if select(plan_path, recipe_path, recipe_sha) != plan:
                    raise ValueError("Plan changed before task")
                result = execute_task(plan, task, recipe_path, recipe_sha, job)
                row["wrapper_status"] = result["status"]
                row["verification"] = record(Path(task["run"]["measurement_directory"]).parent / "verification.json")
                if result["status"] != "command_exited_zero":
                    raise ValueError("Native preparation, measurement or runtime verification failed")
                native = result["measurement"]["native"]
                if native["exit_code"] != 0 or native["timed_out"] is not False:
                    raise ValueError("Native command failed or timed out")
                if select(plan_path, recipe_path, recipe_sha) != plan:
                    raise ValueError("Plan changed after task")
                row.update(status="measurement_completed", native_wall_s=result["measurement"]["native_wall_s"])
            except Exception as error:
                failed = index
                row.update(status="failed", error_type=type(error).__name__, error=str(error))
        save(output / f"task_{index:02d}.json", row)
        rows.append(row)
        save(output / f"progress_{index:02d}.json", dict(tasks=rows, failed_index=failed, job_id=job))
    result = dict(status="root_context_native_completed" if failed is None else "root_context_native_stopped",
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
    raise SystemExit(0 if result["status"] == "root_context_native_completed" else 1)

"""Launch one pinned task of the complete dual-collector overhead experiment."""

import argparse
from functools import partial
import hashlib
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_tools.launch_dgx_native_run import read_pinned, native_enumerator
from benchmark_tools.measure_native_scaling_run import measure_run
from benchmark_tools.measure_native_dual_bracket_step import measure as periodic
from benchmark_tools.measure_frontier_boundary_step import measure as boundary
from benchmark_tools.probe_dgx_step_separation import save

PLAN_SHA = "1160a669a4033c66e7bdde5baddf429e83b9328d3821f99291fd29210cac8ab9"
PROTOCOL_SHA = "7df6f18476d20d9e96ab1a2f1fd69956e42347bf519ba93f35c0edeb1c4e761a"
PURPOSE = "dual_collector_incremental_overhead"


def authorization(recipe_sha):
    return dict(purpose=PURPOSE, execution_authorized=True, scientific_execution_authorized=False,
                plan_sha256=PLAN_SHA, recipe_sha256=recipe_sha, allowed_indices=list(range(18)))


def select(plan_path, recipe_path, recipe_sha, auth_path, auth_sha, index):
    if type(index) is not int or index not in range(18):
        raise ValueError("Require assigned integer task index0-17")
    plan = read_pinned(plan_path, PLAN_SHA)
    auth = read_pinned(auth_path, auth_sha)
    if json.dumps(auth, sort_keys=True) != json.dumps(authorization(recipe_sha), sort_keys=True):
        raise ValueError("Authorization differs from complete engineering scope")
    if ([row["index"] for row in plan["runs"]] != list(range(18)) or plan["purpose"] != PURPOSE
            or plan["protocol_sha256"] != PROTOCOL_SHA or plan["execution_authorized"] is not False
            or plan["scientific_timings_admitted"] is not False):
        raise ValueError("Unexpected complete overhead plan")
    base = Path(__file__).resolve().parent
    protocol = base / "results/DUAL_COLLECTOR_OVERHEAD_PROTOCOL_20260919.md"
    if hashlib.sha256(protocol.read_bytes()).hexdigest() != PROTOCOL_SHA:
        raise ValueError("Protocol changed")
    recipe = read_pinned(recipe_path, recipe_sha)
    files = {row["path"]: row for row in recipe["records"] if row["kind"] == "file"}
    for path in [*base.glob("*.py"), protocol, plan_path.resolve()]:
        if files.get(str(path), {}).get("sha256") != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError("Required source/plan missing from pinned recipe")
    return plan, auth, plan["runs"][index]


def launch(plan_path, recipe_path, recipe_sha, auth_path, auth_sha, index):
    plan, auth, task = select(plan_path, recipe_path, recipe_sha, auth_path, auth_sha, index)
    if (Path(sys.executable).absolute() != Path(plan["launcher_python"])
            or os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "20"
            or os.environ.get("SLURM_MEM_PER_NODE") != "98304"
            or os.environ.get("SLURM_ARRAY_TASK_ID") != str(index)):
        raise ValueError("Require pinned interpreter and matching DGX task/resources")
    cache = str(Path(plan["cache_directory"]) / f"cache_{index}")
    if (not sys.dont_write_bytecode or sys.pycache_prefix != cache or Path(cache).exists()
            or os.environ.get("PYTHONHASHSEED") != "0"):
        raise ValueError("Require fresh isolated cache and frozen Python startup")
    if any(os.environ.get(key) for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT", "PYTHONPATH")):
        raise ValueError("Unexpected loader/Python override")
    if Path("/etc/ld.so.preload").exists():
        raise ValueError("Unexpected system preload")
    run = task["run"]
    for key in plan["unset_environment"]:
        os.environ.pop(key, None)
    os.environ.update(plan["environment_overrides"])
    os.environ["PATH"] = os.pathsep.join(plan["environment_paths"][run["environment_role"]])
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    os.environ["PYTHONPYCACHEPREFIX"] = cache
    runtime = [(item["path"], item["sha256"]) for item in plan["runtime_manifests"]]
    runtime.append((str(recipe_path), recipe_sha))
    collector = periodic if task["mode"] == "periodic" else partial(boundary, native_pressure=True)
    if task["mode"] not in {"periodic", "boundary"}:
        raise ValueError("Unknown collector arm")
    os.chdir(run["cwd"])
    result = measure_run(run, plan["order"], runtime, native_enumerator(plan["enumerator"]),
                         collector, int(os.environ["SLURM_JOB_ID"]), timeout_s=plan["native_timeout_s"])
    if read_pinned(auth_path, auth_sha) != auth or read_pinned(plan_path, PLAN_SHA) != plan:
        raise ValueError("Launch authorization or plan changed")
    save(Path(run["measurement_directory"]).parent / "overhead_task.json",
         dict(task=task, authorization=auth, authorization_sha256=auth_sha,
              recipe_sha256=recipe_sha, plan_sha256=PLAN_SHA, scientific_timings_admitted=False))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "recipe", "authorization"):
        parser.add_argument("--"+name, type=Path, required=True)
    for name in ("recipe-sha", "authorization-sha"):
        parser.add_argument("--"+name, required=True)
    parser.add_argument("--index", type=int, required=True)
    args = parser.parse_args()
    result = launch(args.plan.resolve(), args.recipe.resolve(), args.recipe_sha,
                    args.authorization.resolve(), args.authorization_sha, args.index)
    raise SystemExit(0 if result["status"] == "command_exited_zero" else 1)

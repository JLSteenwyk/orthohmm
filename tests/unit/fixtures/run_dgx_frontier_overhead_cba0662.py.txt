"""Execute one explicitly authorized task from the frozen native overhead panel."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_tools.launch_dgx_native_run import read_pinned, native_enumerator
from benchmark_tools.measure_native_scaling_run import measure_run
from benchmark_tools.measure_native_frontier_step import measure as periodic
from benchmark_tools.measure_frontier_boundary_step import measure as boundary

PLAN_SHA = "58e482e6f123bc82de5e98df3c93f4bf0ca3a47d2302ffcaa47290c7d1e20685"
ROOT = Path("/home/jlsteenwyk/projects/orthohmm-publication")


def select(plan_path, authorization_path, authorization_sha, recipe, recipe_sha, index):
    plan = read_pinned(plan_path, PLAN_SHA)
    auth = read_pinned(authorization_path, authorization_sha)
    if (type(index) is not int or not 0 <= index < 18
            or [row["index"] for row in plan["runs"]] != list(range(18))):
        raise ValueError("Require complete 18-task panel and valid integer index")
    expected = dict(purpose="native_frontier_incremental_overhead", execution_authorized=True,
                    scientific_execution_authorized=False, plan_sha256=PLAN_SHA,
                    recipe_sha256=recipe_sha, allowed_indices=list(range(18)))
    if (auth != expected or type(auth.get("execution_authorized")) is not bool
            or type(auth.get("scientific_execution_authorized")) is not bool
            or any(type(value) is not int for value in auth.get("allowed_indices", []))):
        raise ValueError("Authorization does not match exact engineering scope/recipe")
    if plan["execution_authorized"] is not False or plan["scientific_timings_admitted"] is not False:
        raise ValueError("Unexpected scientific or embedded execution authorization")
    manifest = read_pinned(recipe, recipe_sha)
    required = [Path(__file__).resolve(), Path(periodic.__code__.co_filename).resolve(),
                Path(boundary.__code__.co_filename).resolve(), Path(plan_path).resolve()]
    files = {row["path"]: row for row in manifest["records"] if row["kind"] == "file"}
    for path in required:
        if str(path) not in files or files[str(path)]["sha256"] != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError("Launcher, collector or plan missing from pinned recipe")
    return plan, auth, plan["runs"][index]


def launch(plan_path, authorization_path, authorization_sha, recipe, recipe_sha, index):
    plan, auth, row = select(plan_path, authorization_path, authorization_sha, recipe, recipe_sha, index)
    if (os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "20"
            or os.environ.get("SLURM_MEM_PER_NODE") != "98304"):
        raise ValueError("Require spark-7ff0/20CPU/96GiB")
    if not sys.dont_write_bytecode or os.environ.get("PYTHONHASHSEED") != "0":
        raise ValueError("Require disabled bytecode writes and frozen hash seed")
    if any(os.environ.get(key) for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT")) or Path("/etc/ld.so.preload").exists():
        raise ValueError("Unexpected loader overrides")
    run = row["run"]
    for key in plan["unset_environment"]:
        os.environ.pop(key, None)
    os.environ.update(plan["environment_overrides"])
    os.environ["PATH"] = os.pathsep.join(plan["environment_paths"][run["environment_role"]])
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    os.environ["PYTHONPYCACHEPREFIX"] = str(ROOT / "frontier_overhead_v1" / f"cache_{index}")
    if Path(os.environ["PYTHONPYCACHEPREFIX"]).exists():
        raise ValueError("Cache path exists")
    runtime = [(item["path"], item["sha256"]) for item in plan["runtime_manifests"]]
    runtime.append((str(recipe), recipe_sha))
    collector = {"periodic": periodic, "boundary": boundary}[row["mode"]]
    os.chdir(run["cwd"])
    result = measure_run(run, plan["order"], runtime, native_enumerator(plan["enumerator"]),
                         collector, int(os.environ["SLURM_JOB_ID"]), timeout_s=plan["native_timeout_s"])
    if read_pinned(authorization_path, authorization_sha) != auth:
        raise ValueError("Authorization changed during execution")
    with (Path(run["measurement_directory"]).parent / "overhead_task.json").open("x") as stream:
        json.dump(dict(task=row, authorization=auth, authorization_sha256=authorization_sha,
                       plan_sha256=PLAN_SHA, scientific_timings_admitted=False), stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "authorization", "recipe"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("authorization-sha", "recipe-sha"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--index", type=int, required=True)
    args = parser.parse_args()
    result = launch(args.plan.resolve(), args.authorization.resolve(), args.authorization_sha,
                    args.recipe.resolve(), args.recipe_sha, args.index)
    raise SystemExit(0 if result["status"] == "command_exited_zero" else 1)

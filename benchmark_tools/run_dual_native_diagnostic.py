"""Derive and execute only the three frozen native dual-bracket diagnostics."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_tools.launch_dgx_native_run import read_pinned, native_enumerator
from benchmark_tools.prepare_frontier_overhead_panel import relocate, METHODS
from benchmark_tools.measure_native_scaling_run import measure_run
from benchmark_tools.measure_native_dual_bracket_step import measure

ROOT = Path("/home/jlsteenwyk/projects/orthohmm-publication")
BASE_SHA = "b644e165dbf4d0beabf1cf4d9b6c314de522e3ebd1b91598ebebea99094c8fff"
PROTOCOL_SHA = "488125e6e349a83965a6f8ef75dcc477ce4b8bba5f5a13c0dedb09e88c4e9460"
INDICES = (1, 3, 8)


def build(base_path, protocol_path):
    base = read_pinned(base_path, BASE_SHA)
    if hashlib.sha256(protocol_path.read_bytes()).hexdigest() != PROTOCOL_SHA:
        raise ValueError("Native diagnostic protocol changed")
    runs = []
    for index, (original, method) in enumerate(zip(INDICES, METHODS)):
        template = base["runs"][original]
        if (template["index"] != original or template["method"] != method
                or template["mode"] != "periodic" or template["run"]["proteomes"] != 4
                or template["run"]["dataset"]["proteins"] != 73266):
            raise ValueError("Unexpected frozen native template")
        run = relocate(template["run"], str(ROOT / "pressure_frontier_overhead_v2" / f"run_{original:02d}"),
                       str(ROOT / "dual_bracket_native_v1" / f"run_{index:02d}"))
        run.update(index=index, repeat=0)
        runs.append(dict(index=index, original_index=original, method=method, run=run))
    plan = {k: base[k] for k in ("runtime_manifests", "enumerator", "environment_paths", "environment_overrides",
                                 "unset_environment", "order", "launcher_python", "core_commit")}
    plan.update(status="prospective_dual_bracket_native_diagnostic", runs=runs,
                baseline_sha256=BASE_SHA, protocol_sha256=PROTOCOL_SHA,
                native_timeout_s=900, interval_s=1., execution_authorized=False,
                scientific_timings_admitted=False, publication_ready=False)
    return plan


def select(plan_path, plan_sha, recipe_path, recipe_sha, index):
    if type(index) is not int or index not in range(3):
        raise ValueError("Only three diagnostic indices are permitted")
    folder = Path(__file__).resolve().parent / "results"
    plan = read_pinned(plan_path, plan_sha)
    expected = build(folder / "dgx_pressure_overhead_plan_v2_20260919.json",
                     folder / "DUAL_BRACKET_NATIVE_PROTOCOL_20260919.md")
    if plan != expected:
        raise ValueError("Plan differs from frozen diagnostic derivation")
    recipe = read_pinned(recipe_path, recipe_sha)
    files = {r["path"]: r for r in recipe["records"] if r["kind"] == "file"}
    required = [Path(__file__).resolve(), Path(measure.__code__.co_filename).resolve(), plan_path.resolve(),
                folder / "DUAL_BRACKET_NATIVE_PROTOCOL_20260919.md",
                folder / "dgx_pressure_overhead_plan_v2_20260919.json"]
    for path in required:
        if str(path) not in files or files[str(path)]["sha256"] != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError("Launcher, collector or input missing from pinned recipe")
    return plan, plan["runs"][index]


def launch(plan_path, plan_sha, recipe_path, recipe_sha, index):
    plan, task = select(plan_path, plan_sha, recipe_path, recipe_sha, index)
    if (Path(sys.executable).absolute() != Path(plan["launcher_python"])
            or os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "20"
            or os.environ.get("SLURM_MEM_PER_NODE") != "98304"):
        raise ValueError("Require pinned interpreter and DGX20CPU/96GiB allocation")
    if not sys.dont_write_bytecode or os.environ.get("PYTHONHASHSEED") != "0":
        raise ValueError("Require disabled bytecode and frozen hash seed")
    if any(os.environ.get(k) for k in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT")) or Path("/etc/ld.so.preload").exists():
        raise ValueError("Unexpected loader overrides")
    for key in plan["unset_environment"]:
        os.environ.pop(key, None)
    os.environ.update(plan["environment_overrides"])
    run = task["run"]
    os.environ["PATH"] = os.pathsep.join(plan["environment_paths"][run["environment_role"]])
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    os.environ["PYTHONPYCACHEPREFIX"] = str(ROOT / "dual_bracket_native_v1" / f"cache_{index}")
    if Path(os.environ["PYTHONPYCACHEPREFIX"]).exists():
        raise ValueError("Require fresh cache path")
    runtime = [(r["path"], r["sha256"]) for r in plan["runtime_manifests"]]
    runtime.append((str(recipe_path), recipe_sha))
    os.chdir(run["cwd"])
    result = measure_run(run, plan["order"], runtime, native_enumerator(plan["enumerator"]),
                         measure, int(os.environ["SLURM_JOB_ID"]), timeout_s=900)
    if read_pinned(plan_path, plan_sha) != plan:
        raise ValueError("Plan changed during execution")
    with (Path(run["measurement_directory"]).parent / "dual_diagnostic_task.json").open("x") as stream:
        json.dump(dict(task=task, plan_sha256=plan_sha, recipe_sha256=recipe_sha,
                       purpose="three_native_dual_bracket_diagnostics", scientific_timings_admitted=False),
                  stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--plan-sha")
    parser.add_argument("--recipe", type=Path)
    parser.add_argument("--recipe-sha")
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    if args.prepare:
        folder = Path(__file__).resolve().parent / "results"
        plan = build(folder / "dgx_pressure_overhead_plan_v2_20260919.json",
                     folder / "DUAL_BRACKET_NATIVE_PROTOCOL_20260919.md")
        with args.plan.open("x") as stream:
            json.dump(plan, stream, indent=2, sort_keys=True)
            stream.write("\n")
    else:
        if not args.plan_sha or not args.recipe or not args.recipe_sha or args.index is None:
            parser.error("Execution requires pinned plan, recipe and task index")
        result = launch(args.plan.resolve(), args.plan_sha, args.recipe.resolve(), args.recipe_sha, args.index)
        raise SystemExit(0 if result["status"] == "command_exited_zero" else 1)

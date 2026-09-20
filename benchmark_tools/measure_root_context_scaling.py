"""Per-run measurement composition for the frozen replacement scaling plan.

This library does not authorize or submit runs. Its caller must first verify
execution authorization, environmental policy, allocation and recipe identity.
There is deliberately no standalone execution CLI.
"""

import os
from pathlib import Path
import re

from benchmark_tools.launch_dgx_native_run import read_pinned, native_enumerator
from benchmark_tools.measure_native_root_context import measure_native_run
from benchmark_tools.measure_native_scaling_run import measure_run
from benchmark_tools.prepare_root_context_scaling import OUTPUT_ROOT

PLAN_SHA = "65e0f850f32d09e0049e7e55c8700637588d3b7857aa59a2a70eb942d6c5b348"


def load_task(plan_path, index):
    if type(index) is not int or not 0 <= index < 27:
        raise ValueError("Require original task index in [0, 26]")
    plan = read_pinned(plan_path, PLAN_SHA)
    task = plan["runs"][index]
    if type(task["index"]) is not int or task["index"] != index:
        raise ValueError("Task identity differs")
    matching = [order for order in plan["orders"]
                if order["input_directory"] == task["run"]["dataset"]["input_directory"]]
    if len(matching) != 1:
        raise ValueError("Require one native input order for this dataset")
    return plan, task, matching[0]


def measure_task(plan_path, index, recipe_path, recipe_sha, job):
    if type(job) is not int or job <= 0:
        raise ValueError("Require positive scheduler job identity")
    if not isinstance(recipe_sha, str) or not re.fullmatch(r"[0-9a-f]{64}", recipe_sha):
        raise ValueError("Require recipe SHA-256")
    plan, task, order = load_task(plan_path, index)
    recipe_path = Path(recipe_path).resolve()
    run = task["run"]
    cache = OUTPUT_ROOT / f"cache_{index:02d}"
    if cache.exists() or cache.is_symlink():
        raise ValueError("Task cache prefix must be absent, including dangling links")
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
        return measure_run(run, order, specifications, native_enumerator(plan["enumerator"]),
            measure_native_run, job, cpus=plan["allocation"]["cpus"],
            memory_gib=plan["allocation"]["memory_gib"], timeout_s=plan["native_timeout_s"])
    finally:
        os.chdir(cwd)
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value

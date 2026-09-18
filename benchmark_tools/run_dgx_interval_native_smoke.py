"""Frozen three-method native smoke with complete-command interval observation."""

import argparse
import copy
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from launch_dgx_native_run import read_pinned, native_enumerator
from measure_native_scaling_run import measure_run
from measure_native_interval_step import measure

SPEC_SHA = "8af8a049c69b7e479e8443e1e0e14b74c78d67b345c03d7f557768488d057fa2"
ROOT = Path("/home/jlsteenwyk/projects/orthohmm-publication")


def relocate(value):
    old, new = str(ROOT / "launcher_smoke_v1"), str(ROOT / "interval_native_smoke_v1")
    if isinstance(value, str):
        return new + value[len(old):] if value.startswith(old + "/") else value
    if isinstance(value, dict):
        return {key: relocate(item) for key, item in value.items()}
    if isinstance(value, list):
        return [relocate(item) for item in value]
    return value


def launch(spec_path, index, recipe, recipe_sha):
    if index not in (0, 1, 2):
        raise ValueError("Only three frozen smoke methods")
    if os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "20" or os.environ.get("SLURM_MEM_PER_NODE") != "98304":
        raise ValueError("Require spark-7ff0/20CPU/96GiB")
    if not sys.dont_write_bytecode or os.environ.get("PYTHONHASHSEED") != "0":
        raise ValueError("Require no bytecode writes and frozen hash seed")
    if any(os.environ.get(k) for k in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT")) or Path("/etc/ld.so.preload").exists():
        raise ValueError("Unexpected loader overrides")
    spec = read_pinned(spec_path, SPEC_SHA)
    run = relocate(copy.deepcopy(spec["runs"][index]))
    for key in spec["unset_environment"]:
        os.environ.pop(key, None)
    os.environ.update(spec["environment_overrides"])
    os.environ["PATH"] = os.pathsep.join(spec["environment_paths"][run["environment_role"]])
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    os.environ["PYTHONPYCACHEPREFIX"] = str(ROOT / "interval_native_smoke_v1" / f"cache_{index}")
    if Path(os.environ["PYTHONPYCACHEPREFIX"]).exists():
        raise ValueError("Cache path exists")
    runtime = [(row["path"], row["sha256"]) for row in spec["runtime_manifests"]]
    runtime.append((str(recipe), recipe_sha))
    os.chdir(run["cwd"])
    return measure_run(run, spec["orders"][0], runtime, native_enumerator(spec["enumerator"]),
                       measure, int(os.environ["SLURM_JOB_ID"]), timeout_s=900)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--recipe-sha", required=True)
    parser.add_argument("--index", type=int, required=True)
    args = parser.parse_args()
    result = launch(args.spec, args.index, args.recipe, args.recipe_sha)
    raise SystemExit(0 if result["status"] == "command_exited_zero" else 1)

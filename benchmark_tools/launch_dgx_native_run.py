"""Launch one explicitly authorized, checksum-pinned DGX native run."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import types

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_tools.measure_native_scaling_run import measure_run
from run_verified_slurm_measurement import load_collector


def read_pinned(path, sha):
    path = Path(path)
    contents = path.read_bytes()
    if hashlib.sha256(contents).hexdigest() != sha:
        raise ValueError("Pinned JSON digest differs: " + str(path))
    return json.loads(contents)


def select(spec, index):
    if spec.get("execution_authorized") is not True:
        raise ValueError("Execution has not been explicitly authorized")
    purpose = spec.get("purpose")
    expected = {"launcher_smoke": 3, "scientific_scaling": 27}.get(purpose)
    if expected is None or len(spec["runs"]) != expected:
        raise ValueError("Wrong purpose or complete run inventory")
    if [row["index"] for row in spec["runs"]] != list(range(expected)) or not 0 <= index < expected:
        raise ValueError("Invalid run index/order")
    if purpose == "scientific_scaling":
        proof = spec["validated_launcher_smoke"]
        smoke = read_pinned(proof["path"], proof["sha256"])
        if smoke["status"] != "all_native_launcher_smokes_admitted":
            raise ValueError("Scientific launch requires admitted native launcher smokes")
        baseline = read_pinned(spec["original_plan"]["path"], spec["original_plan"]["sha256"])
        if spec["runs"] != baseline["runs"]:
            raise ValueError("Scientific runs differ from the frozen DGX command plan")
    overhead = read_pinned(spec["overhead"]["path"], spec["overhead"]["sha256"])
    if overhead["status"] != "verified_overhead_panel_evaluated" or overhead["all_protocol_gates_met"] is not True:
        raise ValueError("Verified overhead gate has not passed")
    return spec["runs"][index]


def native_enumerator(source):
    path = Path(source["path"])
    def enumerate_files(directory):
        contents = path.read_bytes()
        if hashlib.sha256(contents).hexdigest() != source["sha256"]:
            raise ValueError("Frozen enumerator source changed")
        module = types.ModuleType("frozen_native_files")
        module.__file__ = str(path)
        exec(compile(contents, str(path), "exec"), module.__dict__)
        return module.fetch_fasta_files(directory)
    return enumerate_files


def launch(spec, index, recipe_manifest, recipe_sha):
    run = select(spec, index)
    if platform.node() != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "20":
        raise ValueError("Require the dedicated20CPU DGX task")
    if spec["resource_plan"]["cpus"] != 20 or spec["resource_plan"]["memory_gib"] != 96:
        raise ValueError("Unexpected matched resource plan")
    cache = Path(run["measurement_directory"]).parent.with_name(Path(run["measurement_directory"]).parent.name + "_cache")
    if (not sys.dont_write_bytecode or sys.pycache_prefix != str(cache) or cache.exists()
            or os.environ.get("PYTHONHASHSEED") != "0"):
        raise ValueError("Require frozen hash seed and absent isolated bytecode cache at interpreter startup")
    for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        if os.environ.get(key):
            raise ValueError("Unexpected dynamic loader override")
    if Path("/etc/ld.so.preload").exists():
        raise ValueError("Unreviewed system preload configuration")
    for key in spec["unset_environment"]:
        os.environ.pop(key, None)
    os.environ.update(spec["environment_overrides"])
    os.environ["PATH"] = os.pathsep.join(spec["environment_paths"][run["environment_role"]])
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    os.environ["PYTHONPYCACHEPREFIX"] = str(cache)
    runtime_specs = [(row["path"], row["sha256"]) for row in spec["runtime_manifests"]]
    runtime_specs.append((str(recipe_manifest), recipe_sha))
    order = [row for row in spec["orders"] if row["input_directory"] == run["dataset"]["input_directory"]]
    if len(order) != 1:
        raise ValueError("Missing or ambiguous frozen native input order")
    os.chdir(run["cwd"])
    def collector(*args, **kwargs):
        measured = load_collector(Path(spec["collector_directory"]), runtime_specs)
        return measured(*args, **kwargs)
    return measure_run(run, order[0], runtime_specs, native_enumerator(spec["enumerator"]), collector,
                       int(os.environ["SLURM_JOB_ID"]), timeout_s=spec["native_timeout_s"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--spec-sha256", required=True)
    parser.add_argument("--recipe-manifest", type=Path, required=True)
    parser.add_argument("--recipe-sha256", required=True)
    parser.add_argument("--index", type=int, required=True)
    args = parser.parse_args()
    spec = read_pinned(args.spec, args.spec_sha256)
    result = launch(spec, args.index, args.recipe_manifest, args.recipe_sha256)
    raise SystemExit(0 if result["status"] == "command_exited_zero" else 1)

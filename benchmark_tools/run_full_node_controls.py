"""Execute the frozen nine-trial engineering panel, never scientific timings."""

import argparse
import hashlib
import os
from pathlib import Path
import platform
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.run_verified_slurm_measurement import run_checked
from benchmark_tools.run_full_node_control_trial import trial
from benchmark_tools.probe_dgx_step_separation import save

ORDER = (("steady", "churn", "contended"), ("churn", "contended", "steady"),
         ("contended", "steady", "churn"))
PROTOCOL_SHA = "76ebeb56b1b7874aa03528064d77fd1c38ec1151a74095204a8393632b0e834d"
PLAN_SHA = "029f19b0e21f356387ae4e2df50310a225cdba4f356037814646af57226e1113"


def summarize(rows):
    expected = [(block, mode) for block, modes in enumerate(ORDER) for mode in modes]
    if [(row["block"], row["mode"]) for row in rows] != expected:
        raise ValueError("Control inventory/order differs")
    valid = all(row["status"] == "workload_validated" for row in rows)
    responses = [dict(block=row["block"], detected=row.get("positive_control_detected")
                 if row["status"] == "workload_validated" else None)
                 for row in rows if row["mode"] == "contended"]
    return dict(all_workloads_valid=valid, positive_controls=responses,
                all_positive_controls_detected=all(row["detected"] is True for row in responses),
                scientific_timings_admitted=False)


def panel(output, job):
    output.mkdir(exist_ok=False)
    rows = []
    for block, modes in enumerate(ORDER):
        for mode in modes:
            directory = output / f"trial_{len(rows):02d}"
            try:
                row = dict(trial(directory, mode, job), block=block)
            except Exception as error:
                row = dict(block=block, mode=mode, job_id=job, status="failed",
                           error_type=type(error).__name__, error=str(error),
                           scientific_timings_admitted=False)
            directory.mkdir(exist_ok=True)
            save(directory / "panel_trial.json", row)
            rows.append(row)
    result = dict(status="full_node_controls_completed", trials=rows, job_id=job,
                  summary=summarize(rows), protocol_sha256=PROTOCOL_SHA,
                  scientific_timings_admitted=False, publication_ready=False)
    save(output / "result.json", result)
    return result


def preflight(recipe_path, recipe_sha):
    base = Path(__file__).resolve().parent
    protocol = base / "results/FULL_NODE_CPU_CONTROL_PROTOCOL_20260919.md"
    if hashlib.sha256(protocol.read_bytes()).hexdigest() != PROTOCOL_SHA:
        raise ValueError("Prospective protocol changed")
    plan_path = base / "results/dgx_dual_native_plan_20260919.json"
    plan = read_pinned(plan_path, PLAN_SHA)
    recipe = read_pinned(recipe_path, recipe_sha)
    files = {row["path"]: row for row in recipe["records"] if row["kind"] == "file"}
    for path in [*base.glob("*.py"), protocol, plan_path]:
        if files.get(str(path), {}).get("sha256") != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError("Recipe does not pin every control source/input")
    if (platform.node() != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "20"
            or os.environ.get("SLURM_MEM_PER_NODE") != "98304"
            or Path(sys.executable).absolute() != Path(plan["launcher_python"])
            or len(os.sched_getaffinity(0)) != 20):
        raise ValueError("Require pinned DGX interpreter and20CPU/96GiB allocation")
    if not sys.dont_write_bytecode or os.environ.get("PYTHONHASHSEED") != "0":
        raise ValueError("Require disabled bytecode and frozen hash seed")
    if not sys.pycache_prefix or Path(sys.pycache_prefix).exists():
        raise ValueError("Require an absent isolated bytecode cache prefix")
    if any(os.environ.get(key) for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT", "PYTHONPATH")):
        raise ValueError("Unexpected loader or Python overrides")
    if Path("/etc/ld.so.preload").exists():
        raise ValueError("Unexpected system preload")
    job = int(os.environ["SLURM_JOB_ID"])
    if job <= 0:
        raise ValueError("Invalid scheduler identity")
    return plan, job


def run(output, recipe_path, recipe_sha):
    plan, job = preflight(recipe_path, recipe_sha)
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    os.environ["PYTHONPYCACHEPREFIX"] = sys.pycache_prefix
    specs = [(row["path"], row["sha256"]) for row in plan["runtime_manifests"]]
    specs.append((str(recipe_path), recipe_sha))
    result = run_checked(specs, output, lambda directory: panel(directory, job))
    save(output / "launch.json", dict(job_id=job, recipe_sha256=recipe_sha,
        protocol_sha256=PROTOCOL_SHA, plan_sha256=PLAN_SHA, executable=sys.executable,
        python=sys.version, kernel=platform.release(), host=platform.node()))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--recipe-sha", required=True)
    args = parser.parse_args()
    result = run(args.output.resolve(), args.recipe.resolve(), args.recipe_sha)
    summary = result.get("measurement", {}).get("summary", {})
    raise SystemExit(0 if result["status"] == "full_node_controls_completed"
                     and summary.get("all_workloads_valid") is True
                     and summary.get("all_positive_controls_detected") is True else 1)

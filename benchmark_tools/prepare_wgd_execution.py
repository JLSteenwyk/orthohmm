"""Bind the frozen WGD command panel to explicit runtime and launcher snapshots."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.run_wgd_application import pinned
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools.snapshot_runtime_trees import inventory


def prepare(repo, launcher, recipe):
    results = repo / "benchmark_tools/results"
    work = repo / "benchmarks/work"
    plan_record = record(results / "biological_wgd_commands_20260917.json")
    if plan_record["sha256"] != "9a4f1af732d7e535818afe8d0da93743da99401dfe37f97433d94a860eaed1f3":
        raise ValueError("Changed native command plan")
    plan = pinned(plan_record)
    subprocess.run(["git", "-C", str(launcher), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    commit = subprocess.check_output(["git", "-C", str(launcher), "rev-parse", "HEAD"], text=True).strip()
    roots = [launcher / "benchmark_tools" / name for name in
             ("run_wgd_application.py", "run_wgd_application.slurm", "snapshot_runtime_trees.py",
              "snapshot_orthohmm_input_order.py", "__init__.py")]
    recipe_data = inventory(roots)
    with recipe.open("x") as handle:
        json.dump(recipe_data, handle, indent=2, sort_keys=True)
        handle.write("\n")
    runtime = record(work / "biological_wgd_runtime_trees_v1.json")
    system = record(work / "biological_wgd_system_trees_v1.json")
    if runtime["sha256"] != "70eda6198d28d6c36c697d2862911a0853d9482c9e15f504f17574ca73f99071":
        raise ValueError("Changed original runtime inventory")
    if system["sha256"] != "5a19ac94e14aff479c7bd7a09a757a8471056d23a9b51fbbc8c51a2d8aee0e27":
        raise ValueError("Changed system runtime inventory")
    software = repo.parents[1] / "SOFTWARE"
    paths = [software / "mafft-7.525-with-extensions/bin", software / "FastTree_v220",
             software / "diamond-linux64", Path("/home/bizon/anaconda3/bin"), Path("/usr/bin"), Path("/bin")]
    return {"purpose": "wgd_application", "execution_authorized": True,
            "command_plan": plan_record, "launcher_commit": commit, "launcher_root": str(launcher),
            "source": record(__file__), "runtime_manifests": [runtime, system, record(recipe)],
            "input_order": record(results / "biological_wgd_input_order_20260917.json"),
            "timeout_s": 85800,
            "environment": {"PATH": ":".join(map(str, paths)), "PYTHONPATH": plan["core_root"],
                            "PYTHONHASHSEED": "0", "PYTHONUNBUFFERED": "1",
                            "PYTHONDONTWRITEBYTECODE": "1", "OMP_NUM_THREADS": "1",
                            "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                            "PYTHONNOUSERSITE": "1"},
            "unset_environment": ["CONDA_PREFIX", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT",
                                  "GOMP_CPU_AFFINITY", "OMP_PROC_BIND", "OMP_PLACES", "OMP_DYNAMIC"],
            "limitations": ["Execution authorization does not admit native outputs or scientific claims.",
                            "Original shared host: no controlled cross-method timing comparison.",
                            "Explicit runtime identities are not a hermetic operating-system snapshot.",
                            "System Qt default.conf broken link retained; no Qt workflow is used.",
                            "Before/after identity cannot exclude temporary changes during inference."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "launcher", "recipe", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.recipe.exists():
        raise FileExistsError("Refusing to replace execution artifacts")
    result = prepare(args.repo.resolve(), args.launcher.resolve(), args.recipe.resolve())
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")

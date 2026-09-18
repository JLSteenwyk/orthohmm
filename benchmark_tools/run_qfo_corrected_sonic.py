"""Run frozen corrected-input SonicParanoid without admitting accuracy outputs."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.prepare_qfo_corrected_sonic import command, validate_runtime, PRIMARY_SHA
from benchmark_tools.snapshot_runtime_trees import verify as verify_tree

PLAN_SHA = "188efb1f18e8f0bd90967670d18015ad500a9c47480333b0b380dff7c10c4c59"
RUNTIME_SHA = "b4ed57defaba417f6627b671df958a2d533f4438bb8ff35a0b95cb820fd61bc9"
INVENTORIES = (
    ("biological_wgd_runtime_trees_v1.json", "70eda6198d28d6c36c697d2862911a0853d9482c9e15f504f17574ca73f99071"),
    ("biological_wgd_system_trees_v1.json", "5a19ac94e14aff479c7bd7a09a757a8471056d23a9b51fbbc8c51a2d8aee0e27"),
)


def environment(plan, runtime):
    root = Path(plan["output_root"])
    return {"HOME": str(Path.home()), "USER": "bizon", "LOGNAME": "bizon",
            "PATH": runtime["path"], "LANG": "C", "LC_ALL": "C",
            "PYTHONHASHSEED": "0", "PYTHONDONTWRITEBYTECODE": "1", "PYTHONNOUSERSITE": "1",
            "PYTHONPYCACHEPREFIX": str(root.with_name(root.name + "_pycache")),
            "PYTHONUNBUFFERED": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}


def compare_resolution(expected, observed):
    validate_runtime(observed)
    for key in ("version", "default_mode", "python", "package_inventory"):
        if observed[key] != expected[key]:
            raise ValueError("SonicParanoid runtime drift: " + key)
    for name, tool in expected["tools"].items():
        if observed["tools"][name]["file"] != tool["file"]:
            raise ValueError("SonicParanoid dependency resolution drift: " + name)


def verify_corrected_inputs(primary, inputs):
    for item in primary["inputs"]:
        check(item)
    directory = Path(primary["input_directory"])
    stage = json.loads((directory / "staging_manifest.json").read_text())
    if len(inputs) != 78 or inputs != stage["input_fastas"]:
        raise ValueError("Changed corrected input inventory")
    expected = {Path(item["path"]).name for item in inputs} | {"staging_manifest.json"}
    if {p.name for p in directory.iterdir()} != expected:
        raise ValueError("Changed corrected input directory membership")
    for item in inputs:
        if Path(item["path"]).parent != directory:
            raise ValueError("Changed input location")
        check(item)


def verify(plan_path):
    plan = read_frozen(plan_path, PLAN_SHA)
    if plan["status"] != "corrected_sonic_command_frozen_unrun" or plan["execution_authorized"] is not False:
        raise ValueError("Unexpected command plan state")
    root = Path(plan["output_root"])
    if plan["native_argv"] != command("/home/bizon/anaconda3/bin/sonicparanoid", root):
        raise ValueError("Changed native command")
    if plan["cwd"] != str(root) or plan["copy_inputs_to"] != str(root / "input"):
        raise ValueError("Changed native working/input directory")
    for item in [plan["source"], *plan["checked_records"]]:
        check(item)
    primary = read_frozen(Path(plan["primary_manifest"]["path"]), PRIMARY_SHA)
    verify_corrected_inputs(primary, plan["input_fastas"])
    runtime = read_frozen(Path(plan["runtime_manifest"]["path"]), RUNTIME_SHA)
    validate_runtime(runtime)
    repo = plan_path.resolve().parents[2]
    inventories = []
    for name, digest in INVENTORIES:
        path = repo / "benchmarks/work" / name
        inventories.append({"manifest": record(path), "verification": verify_tree(read_frozen(path, digest))})
    if Path("/etc/ld.so.preload").exists():
        raise ValueError("Unreviewed loader preload")
    env = environment(plan, runtime)
    code = ("import importlib.util,json,sys; "
            "s=importlib.util.spec_from_file_location('probe',sys.argv[1]); "
            "m=importlib.util.module_from_spec(s); s.loader.exec_module(m); print(json.dumps(m.inspect()))")
    done = subprocess.run(["/home/bizon/anaconda3/bin/python", "-B", "-s", "-c", code, runtime["source"]["path"]],
                          env=env, cwd=Path(plan["input_directory"]).parent,
                          capture_output=True, text=True, check=True, timeout=240)
    compare_resolution(runtime, json.loads(done.stdout))
    if Path(env["PYTHONPYCACHEPREFIX"]).exists():
        raise ValueError("Isolated bytecode prefix must remain absent")
    return plan, env, inventories


def run(plan_path, check_only=False):
    if not check_only and (os.environ.get("SLURM_CPUS_PER_TASK") != "32" or not os.environ.get("SLURM_JOB_ID")):
        raise ValueError("Require scheduled 32-CPU allocation")
    plan, env, inventories = verify(plan_path)
    root = Path(plan["output_root"])
    if root.exists():
        raise FileExistsError("Existing output; no implicit restart")
    if check_only:
        return {"status": "preflight_passed_no_inference", "inventories": inventories}
    root.mkdir(parents=True, exist_ok=False)
    status = root / "execution.json"
    report = {"status": "preparing", "source": record(__file__), "plan": record(plan_path),
              "runtime_before": inventories, "job_id": os.environ["SLURM_JOB_ID"],
              "node": os.uname().nodename, "native_argv": plan["native_argv"], "cwd": str(root),
              "environment": env, "accuracy_admitted": False, "native_outputs_validated": False,
              "limitation": "Shared-host inference, not dedicated timing; native output and scoring admission remain separate."}
    try:
        directory = root / "input"
        directory.mkdir(exist_ok=False)
        copies = []
        for item in plan["input_fastas"]:
            check(item)
            target = directory / Path(item["path"]).name
            shutil.copy2(item["path"], target)
            copied = record(target)
            if any(copied[k] != item[k] for k in ("sha256", "bytes")):
                raise ValueError("Input copy differs")
            copies.append(copied)
        report.update(status="running", copied_inputs=copies, started_epoch=time.time())
        status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        with (root / "native.log").open("xb") as log:
            done = subprocess.run(["/usr/bin/time", "-v", "-o", str(root / "time.txt"), *plan["native_argv"]],
                                  cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT)
        report.update(exit_code=done.returncode, finished_epoch=time.time())
        _, _, report["runtime_after"] = verify(plan_path)
        check(report["plan"])
        for item in copies:
            check(item)
        if {p.name for p in directory.glob("*.fasta")} != {Path(r["path"]).name for r in copies}:
            raise ValueError("Changed copied FASTA membership")
        report.update(outputs=[record(p) for p in sorted((root / "output").rglob("*")) if p.is_file()],
                      log=record(root / "native.log"), timing=record(root / "time.txt"))
        if done.returncode:
            raise RuntimeError(f"Native process failed: {done.returncode}")
        report["status"] = "process_succeeded_pending_native_admission"
    except Exception as exc:
        report.update(status="failed", error=str(exc), finished_epoch=time.time())
        raise
    finally:
        status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    result = run(args.plan.resolve(), args.check_only)
    print(json.dumps({"status": result["status"]}))

"""Freeze runtime state and execute fresh corrected-input Proteinortho inference."""

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
from benchmark_tools.prepare_qfo_corrected_proteinortho import native_command, PRIMARY_SHA
from benchmark_tools.run_qfo_corrected_primary import verify as verify_primary

PLAN_SHA = "65f9baedeb6f6cbb66198e2a6a24ae8531b48eea9a04809c2c2443d570d4f93a"
CONFIG = Path("/usr/local/etc/singularity")
HELPERS = Path("/usr/local/libexec/singularity/bin")
PROBE = '''import glob, hashlib, json, os, shutil
paths = set(glob.glob('/usr/local/bin/proteinortho*'))
for name in ('proteinortho', 'diamond', 'perl', 'python3', 'sort', 'rm', 'mkdir'):
    path = shutil.which(name)
    if path is None:
        raise RuntimeError('Missing dependency: ' + name)
    paths.add(path)
files = []
for path in sorted(paths):
    resolved = os.path.realpath(path)
    if os.path.isfile(resolved):
        with open(resolved, 'rb') as stream:
            digest = hashlib.sha256(stream.read()).hexdigest()
        files.append({'path': path, 'resolved': resolved, 'bytes': os.path.getsize(resolved), 'sha256': digest})
print(json.dumps({'files': files, 'environment': {k:v for k,v in os.environ.items()
    if k in ('PATH','LD_LIBRARY_PATH','PERL5LIB','PERL5OPT','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','LANG','LC_ALL')}}))
'''


def native_environment(plan):
    # Do not propagate scheduler-exported container overrides or language startup hooks.
    observed = plan["observed_environment"]
    env = {"HOME": str(Path.home()), "USER": "bizon", "LOGNAME": "bizon",
           "PATH": observed["PATH"], "LANG": "C", "LC_ALL": "C"}
    for key in ("LD_LIBRARY_PATH", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        if key in observed:
            env[key] = observed[key]
    return env


def runtime_snapshot(plan, env, cwd):
    paths = sorted(p for root in (CONFIG, HELPERS) for p in root.rglob("*") if p.is_file())
    files = [record(p) for p in paths]
    command = plan["native_argv"][:5] + ["python3", "-c", PROBE]
    done = subprocess.run(command, env=env, cwd=cwd, text=True, capture_output=True, check=True, timeout=120)
    for item in files:
        check(item)
    return {"host_files": files, "container": json.loads(done.stdout), "environment": env,
            "limitation": "Effective dependency/configuration snapshot, not a complete system-library inventory or execution trace."}


def verify_plan(path):
    plan = read_frozen(path, PLAN_SHA)
    if plan["status"] != "corrected_proteinortho_command_frozen_unrun" or plan["execution_authorized"] is not False:
        raise ValueError("Unexpected plan state")
    for item in [plan["source"], *plan["checked_records"]]:
        check(item)
    primary = read_frozen(Path(plan["primary_manifest"]["path"]), PRIMARY_SHA)
    _, inputs = verify_primary(primary)
    if inputs != plan["input_fastas"]:
        raise ValueError("Changed corrected input inventory")
    expected = native_command(plan["runtime"]["path"], plan["image"]["path"],
                              [Path(r["path"]).name for r in inputs])
    if expected != plan["native_argv"]:
        raise ValueError("Changed native command")
    return plan


def freeze(plan_path, destination):
    if destination.exists():
        raise FileExistsError(destination)
    plan = verify_plan(plan_path)
    env = native_environment(plan)
    root = Path(plan["input_directory"]).parent
    snapshot = runtime_snapshot(plan, env, root)
    report = {"status": "proteinortho_runtime_frozen_unrun", "source": record(__file__),
              "plan": record(plan_path), "snapshot": snapshot,
              "environment_policy": "Explicit allowlist, fixed C locale; inherited container overrides and language startup hooks excluded.",
              "accuracy_admitted": False}
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


def run(plan_path, runtime_path, runtime_sha):
    if os.environ.get("SLURM_CPUS_PER_TASK") != "32" or not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Require scheduled 32-CPU allocation")
    plan = verify_plan(plan_path)
    frozen = read_frozen(runtime_path, runtime_sha)
    if frozen["status"] != "proteinortho_runtime_frozen_unrun" or frozen["plan"] != record(plan_path):
        raise ValueError("Changed runtime plan binding")
    if record(__file__)["sha256"] != frozen["source"]["sha256"]:
        raise ValueError("Runner differs from frozen runtime source")
    env = native_environment(plan)
    root = Path(plan["output_root"])
    if root.exists():
        raise FileExistsError("Existing output; no implicit restart")
    if runtime_snapshot(plan, env, Path(plan["input_directory"]).parent) != frozen["snapshot"]:
        raise ValueError("Runtime drift before execution")
    root.mkdir(parents=True, exist_ok=False)
    status = root / "execution.json"
    report = {"status": "preparing", "source": record(__file__), "plan": record(plan_path),
              "runtime": record(runtime_path), "job_id": os.environ["SLURM_JOB_ID"],
              "node": os.uname().nodename, "native_argv": plan["native_argv"], "cwd": plan["cwd"],
              "accuracy_admitted": False, "native_outputs_validated": False,
              "limitation": "Shared-host inference is not dedicated timing. Conversion/scoring require separate admission."}
    try:
        directory = Path(plan["cwd"])
        if directory != root / "input":
            raise ValueError("Unexpected native working directory")
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
        if runtime_snapshot(plan, env, directory) != frozen["snapshot"]:
            raise ValueError("Runtime drift in native working directory")
        report.update(status="running", copied_inputs=copies, started_epoch=time.time())
        status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        with (root / "native.log").open("xb") as log:
            done = subprocess.run(["/usr/bin/time", "-v", "-o", str(root / "time.txt"), *plan["native_argv"]],
                                  env=env, cwd=directory, stdout=log, stderr=subprocess.STDOUT)
        report.update(exit_code=done.returncode, finished_epoch=time.time())
        verify_plan(plan_path)
        check(report["runtime"])
        for item in copies:
            check(item)
        if runtime_snapshot(plan, env, directory) != frozen["snapshot"]:
            raise ValueError("Runtime drift after execution")
        report.update(outputs=[record(p) for p in sorted(directory.rglob("*")) if p.is_file()],
                      log=record(root / "native.log"), timing=record(root / "time.txt"))
        if done.returncode:
            raise RuntimeError(f"Native process failed: {done.returncode}")
        report["status"] = "process_succeeded_pending_native_admission"
    except Exception as exc:
        report.update(status="failed", error=str(exc), finished_epoch=time.time())
        raise
    finally:
        status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("freeze", "run"))
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--runtime", type=Path, required=True)
    parser.add_argument("--runtime-sha256")
    args = parser.parse_args()
    if args.mode == "freeze":
        freeze(args.plan.resolve(), args.runtime.resolve())
    else:
        if not args.runtime_sha256:
            parser.error("run requires --runtime-sha256")
        run(args.plan.resolve(), args.runtime.resolve(), args.runtime_sha256)

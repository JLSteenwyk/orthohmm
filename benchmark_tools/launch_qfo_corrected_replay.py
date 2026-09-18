"""Authorize a prepared corrected replay, then exec its unchanged frozen runner."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.checked_replay_payload_worker import corrected_evidence
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_replay import command_for, validate_admission, PLAN_SHA, METHOD
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR = "188fde21860a70da55fee1358177485c970a11a3"
ADMITTER = "7b5214a5d2169338c4cbe80293fae491ee5b951f"
STAGES = ["multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"]


def completed(job):
    text = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    scheduler = require_completed_job(text, job)
    if any(scheduler.get(k) != v for k, v in {"NodeList": "bizon", "AllocCPUS": "2", "ReqMem": "64G"}.items()):
        raise ValueError("Wrong preparation/admission allocation")
    return scheduler, text


def frozen(path, commit):
    if subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip() != commit:
        raise ValueError("Changed frozen replay/preparation/admission executor")
    subprocess.run(["git", "-C", str(path), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)


def validate_plan(root, plan, admission, primary, executor, admitter):
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    output = root / "benchmarks/results/qfo_corrected_checked_replay_v1"
    admission_path = root / "benchmarks/work/qfo_corrected_high_sensitivity_admission_20260918.json"
    primary_path = root / "benchmark_tools/results/qfo_corrected_primary_commands_20260918.json"
    if (plan["source"] != record(executor / "benchmark_tools/prepare_qfo_corrected_replay.py")
            or admission["source"] != record(admitter / "benchmark_tools/admit_qfo_corrected_high_sensitivity.py")
            or plan["admission"] != record(admission_path) or plan["primary_plan"] != record(primary_path)):
        raise ValueError("Changed prepared replay source/admission/primary binding")
    checkpoint, sha = validate_admission(admission, primary, plan["input_fastas"])
    expected = command_for(Path(sys.executable), launcher, output, Path(primary["input_directory"]), checkpoint, sha)
    environment = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                   "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if (plan["native_command"] != expected or primary["methods"][METHOD]["native_argv"][0] != sys.executable
            or plan["output_root"] != str(output) or plan["cwd"] != str(launcher)
            or plan["environment_overrides"] != environment or plan["expected_stages"] != STAGES):
        raise ValueError("Prepared replay differs from frozen scientific settings/paths")
    if output.exists() or output.is_symlink():
        raise FileExistsError("Replay already attempted; no implicit retry")


def prepare(root, preparation_job, admission_job):
    preparation_scheduler, preparation_accounting = completed(preparation_job)
    admission_scheduler, admission_accounting = completed(admission_job)
    executor = root / "benchmarks/work/publication_qfo_corrected_replay_v1"
    admitter = root / "benchmarks/work/publication_qfo_corrected_hmm_admission_v1"
    frozen(executor, EXECUTOR)
    frozen(admitter, ADMITTER)
    plan_path = root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json"
    plan_record = record(plan_path)
    plan, checked_plan, admission_record, _ = corrected_evidence(plan_path, plan_record["sha256"])
    if checked_plan != plan_record:
        raise ValueError("Prepared plan changed during authorization")
    primary_path = root / "benchmark_tools/results/qfo_corrected_primary_commands_20260918.json"
    primary = read_frozen(primary_path, PLAN_SHA)
    admission = json.loads(Path(admission_record["path"]).read_text())
    validate_plan(root, plan, admission, primary, executor, admitter)
    command = [sys.executable, str(executor / "benchmark_tools/run_qfo_corrected_replay.py"),
               "--root", str(root), "--plan", str(plan_path), "--plan-sha256", plan_record["sha256"]]
    checked = [plan_record, admission_record, plan["primary_plan"], plan["source"], admission["source"],
               *plan["checked_records"], *[record(executor / "benchmark_tools" / name) for name in (
                   "run_qfo_corrected_replay.py", "checked_replay_payload_worker.py", "checked_replay_interceptor.py")],
               record(__file__)]
    for item in checked:
        check(item)
    return {"status": "corrected_replay_dispatch_authorized", "source": record(__file__),
            "plan": plan_record, "preparation_scheduler": preparation_scheduler,
            "preparation_accounting": preparation_accounting, "admission_scheduler": admission_scheduler,
            "admission_accounting": admission_accounting, "command": command, "checked_records": checked,
            "execution_authorized": True, "accuracy_evaluated": False, "publication_ready": False,
            "limitations": ["Authorizes only the frozen checked runner, never the unchecked native command.",
                            "Runner must still verify scientific runtime, all four clustering boundaries and output coverage.",
                            "Cached shared-host replay is not matched end-to-end timing; admission and scoring remain required."]}


def launch(root, preparation_job, admission_job, check_only=False):
    result = prepare(root, preparation_job, admission_job)
    if check_only:
        return result
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "32"
            or os.environ.get("SLURM_MEM_PER_NODE") != "196608" or os.uname().nodename != "bizon"):
        raise ValueError("Require scheduled 32-CPU/192-GiB replay on bizon")
    destination = root / "benchmarks/work/qfo_corrected_replay_dispatch_20260918.json"
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    result["job_id"] = os.environ["SLURM_JOB_ID"]
    with destination.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.execv(result["command"][0], result["command"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--preparation-job", type=int, required=True)
    parser.add_argument("--admission-job", type=int, required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    result = launch(args.root.resolve(), args.preparation_job, args.admission_job, args.check_only)
    if args.check_only:
        print(json.dumps({"status": result["status"], "command": result["command"]}))

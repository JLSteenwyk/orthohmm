"""Run only prespecified CPM changes after independently admitted control replay."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

# Spawn replays this file as __mp_main__; preserve the parent's verified
# native import path instead of shadowing it with this orchestration tree.
if __name__ != "__mp_main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

ARMS = ("cpm_low", "cpm_high")
ADMISSION_JOB = "21958"
ADMISSION_COMMIT = "40806f569e4119d3b763fca8d80e5b2ec2685c0e"
ADMISSION_SHA = "fa1252387a5afb3abad905e3bc86ec96bfbf88954d59b1b140fb6a93d1e45762"


def control_evidence(root, control_context):
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.verify_ygob_validation import require_completed_job
    accounting = subprocess.check_output(["sacct", "-j", ADMISSION_JOB, "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, ADMISSION_JOB)
    if scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "2":
        raise ValueError("Wrong control admission allocation")
    executor = root / "benchmarks/work/publication_qfo_cpm_control_admission_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMISSION_COMMIT:
        raise ValueError("Control admission executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    source = record(executor / "benchmark_tools/admit_qfo_cpm_control.py")
    if source["sha256"] != ADMISSION_SHA:
        raise ValueError("Control admission source changed")
    item = record(root / f"benchmarks/work/qfo_cpm_control_admission_{ADMISSION_JOB}.json")
    report = read_frozen(Path(item["path"]), item["sha256"])
    if (report["status"] != "cpm_control_reproduced_and_admitted" or report["source"] != source
            or report["context"] != control_context or report["changed_arms_authorized"] != list(ARMS)
            or report["accuracy_evaluated"] is not False or report["publication_ready"] is not False):
        raise ValueError("Missing or changed control authorization")
    records = [source, item, *report["checked_records"]]
    for entry in records:
        check(entry)
    return {"scheduler": scheduler, "accounting": accounting, "report": report,
            "report_record": item, "checked_records": records, "executor": str(executor)}


def run(root, index):
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.checked_replay_payload_worker import corrected_evidence
    from benchmark_tools.cpm_replay_context import evidence, REPLAY_SHA
    from benchmark_tools.run_qfo_corrected_replay import validate_completion
    from benchmark_tools.audit_historical_profile_ablation import read_partition
    from benchmark_tools.verify_qfo_replay_launcher import verify
    if type(index) is not int or index not in range(len(ARMS)):
        raise ValueError("Unknown prespecified CPM index")
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "32"
            or os.environ.get("SLURM_JOB_NODELIST") != "bizon"
            or os.environ.get("SLURM_ARRAY_TASK_ID") != str(index)):
        raise ValueError("Require scheduled 32-CPU bizon CPM array task")
    arm = ARMS[index]
    plan_path = root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json"
    plan, plan_record, _, names_record = corrected_evidence(plan_path, REPLAY_SHA)
    context = evidence(root, plan, plan_record, arm)
    control_context = evidence(root, plan, plan_record, "control")
    authorization = control_evidence(root, control_context)
    output = Path(context["output_root"])
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    executor = Path(__file__).resolve().parent.parent
    source = record(__file__)
    helpers = [record(p) for p in sorted((executor / "benchmark_tools").glob("*.py"))]
    launcher = Path(context["cwd"])
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    if runtime != plan["runtime"]:
        raise ValueError("Changed frozen CPM runtime")
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "running", "arm": arm, "index": index, "context": context,
        "source": source, "helpers": helpers, "authorization": authorization, "runtime_before": runtime,
        "job_id": os.environ["SLURM_JOB_ID"], "started_epoch": time.time(),
        "executor_commit": subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip(),
        "accuracy_evaluated": False, "publication_ready": False}
    try:
        # Reproduce authorization with the independently frozen validator before inference.
        admission_output = output / "fresh_control_admission.json"
        admission_command = [sys.executable, "-B", str(Path(authorization["executor"]) / "benchmark_tools/admit_qfo_cpm_control.py"),
                             "--root", str(root), "--output", str(admission_output)]
        report["admission_command"] = admission_command
        with (output / "control_validation.log").open("x") as log:
            subprocess.run(admission_command, check=True, stdout=log, stderr=subprocess.STDOUT)
        fresh_record = record(admission_output)
        if read_frozen(admission_output, fresh_record["sha256"]) != authorization["report"]:
            raise ValueError("Fresh control admission disagrees")
        report["fresh_control_admission"] = fresh_record
        command = [sys.executable, "-B", str(Path(__file__).resolve()), "--root", str(root),
                   "--index", str(index), "--worker"]
        report["worker_command"] = command
        with (output / "replay.log").open("x") as log:
            process = subprocess.run(["/usr/bin/time", "-v", "-o", str(output / "time.txt"), *command],
                cwd=launcher, env={**os.environ, **context["environment_overrides"]}, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = process.returncode
        if process.returncode:
            raise RuntimeError("CPM variant failed; no implicit retry")
        worker = json.loads((output / "checked_worker.json").read_text())
        replay = json.loads((output / "replay.json").read_text())
        if worker["context"] != context:
            raise ValueError("Variant worker context changed")
        validate_completion(worker, replay, context["expected_stages"])
        names = Path(names_record["path"]).read_text().splitlines()
        universe = set(names)
        if len(names) != 984137 or len(universe) != len(names):
            raise ValueError("Wrong variant gene universe")
        coverage = []
        for stage in replay["stages"]:
            check(stage["output"])
            groups = read_partition(Path(stage["output"]["path"]), universe)
            coverage.append({"label": stage["label"], "groups": len(groups), "output": stage["output"]})
        report["runtime_after"] = verify(core, launcher, runtime_path)
        if report["runtime_after"] != runtime:
            raise ValueError("CPM runtime changed during replay")
        corrected_evidence(plan_path, REPLAY_SHA)
        evidence(root, plan, plan_record, arm)
        for item in [source, names_record, fresh_record, *helpers, *authorization["checked_records"],
                     *context["checked_records"], *[r["output"] for r in coverage]]:
            check(item)
        report.update(status="cpm_variant_replay_complete_pending_admission", coverage=coverage,
            worker=record(output / "checked_worker.json"), replay=record(output / "replay.json"),
            limitations=["No accuracy evaluation or default promotion; independent replay admission remains required.",
                         "Candidate construction, phylogeny and scoring remain required.",
                         "Cached shared-host runtime is not controlled efficiency evidence."])
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_epoch"] = time.time()
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(2), required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from run_qfo_cpm_control import worker
        worker(args.root.resolve(), ARMS[args.index])
    else:
        run(args.root.resolve(), args.index)

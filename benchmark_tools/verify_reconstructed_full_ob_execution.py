"""Verify execution of job 22376; scientific admission remains a separate step."""

import argparse
import csv
import io
from pathlib import Path
import subprocess

from benchmark_tools.admit_integrated_full_ob import (
    expected_stages, read, verify_private_inputs, verify_stages,
)
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_integrated_full_job import validate
from benchmark_tools.run_integrated_publication_workflow import install_commands
from benchmark_tools.run_publication_pipeline import save
from benchmark_tools.verify_ygob_validation import require_completed_job

JOB = 22376
PLAN_SHA = "cea0dd8efa459a01c57ca0450798005c5de3aeef15a6fda19c181acd35918bee"
RESULTS = Path(__file__).resolve().parent / "results"


def scheduler_gate(accounting, job):
    if job != JOB:
        raise ValueError("Wrong prespecified job")
    row = require_completed_job(accounting, job)
    if (row["AllocCPUS"] != "32" or row["NodeList"] != "bizon"
            or row["ReqMem"] not in {"128G", "128Gn"}):
        raise ValueError("Wrong scheduler resources")
    steps = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    for step in steps:
        if step["JobIDRaw"].startswith(str(job) + "."):
            if step["State"] != "COMPLETED" or step["ExitCode"] != "0:0":
                raise ValueError("Unsuccessful scheduler step")
    return row


def execution_binding(directory, plan, submitted, started, executed):
    pinned = record(directory / "plan.json")
    if (submitted["job_id"] != JOB or submitted["protocol_commit"] != "83f0ebea"
            or submitted["plan"] != pinned
            or started != dict(plan=pinned, job_id=str(JOB))
            or executed["status"] != "integrated_complete_pending_independent_admission"
            or executed["returncode"] != 0 or executed["job_id"] != str(JOB)
            or executed["plan"] != pinned
            or executed["command"] != ["/usr/bin/time", "-v", "-o", str(directory / "time.txt"), *plan["command"]]):
        raise ValueError("Submission or execution binding differs")


def verify(directory, job=JOB):
    if job != JOB:
        raise ValueError("Wrong prespecified job")
    accounting = subprocess.check_output([
        "sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,NodeList,ReqMem"], text=True)
    scheduler = scheduler_gate(accounting, job)
    plan = validate(directory / "plan.json", PLAN_SHA)
    submission = RESULTS / "reconstructed_full_ob_submission_22376.json"
    preparation = RESULTS / "reconstructed_full_ob_preparation_20260929.json"
    submitted = read(submission)
    execution_binding(directory, plan, submitted, read(directory / "execution_started.json"),
                      read(directory / "execution.json"))
    if read(preparation)["plan"] != record(directory / "plan.json"):
        raise ValueError("Preparation plan differs")
    check(submitted["script"])
    executed = read(directory / "execution.json")
    root = directory / "run"
    if any((root / name).exists() for name in ("failure.json", "native/failure.json")):
        raise ValueError("Failed workflow")
    complete = read(root / "complete.json")
    if (executed["complete"] != record(root / "complete.json")
            or complete["status"] != "integrated_install_inference_readback_scoring_complete"
            or complete["dataset"] != "orthobench"
            or complete["native_checkpoint_reuse"] is not False
            or complete["controlled_timing"] is not False
            or complete["started"] != record(root / "started.json")
            or complete["score"] != record(root / "score.json")):
        raise ValueError("Completion binding differs")
    controller = Path(plan["command"][4])
    if record(Path(install_commands.__code__.co_filename))["sha256"] != record(controller)["sha256"]:
        raise ValueError("Controller helper differs")
    stages = expected_stages(plan, directory)
    checked = verify_stages(root, complete, stages)
    started = read(root / "started.json")
    pinned = {r["path"]: r for r in plan["checked_records"]}
    if (started["source"] != record(controller) or started["cpu"] != 32
            or started["dataset"] != "orthobench" or started["attempts"] != 1
            or started["native_checkpoint_reuse"] is not False or not started["inputs"]
            or any(pinned.get(r["path"]) != r for r in started["inputs"])):
        raise ValueError("Workflow start binding differs")
    native = read(root / "native/started.json")
    command = stages[6][1]
    if (native["command"] != command[3:] or native["executable"] != command[0]
            or native["attempts"] != 1 or native["checkpoint_reuse"] is not False):
        raise ValueError("Native execution binding differs")
    data = Path(plan["command"][plan["command"].index("--data") + 1])
    verify_private_inputs(read(data), read(root / "data.json"), root)
    logs = [record(directory / name) for name in ("workflow.log", "time.txt")]
    if executed["logs"] != logs:
        raise ValueError("Execution logs differ")
    checked.extend(logs + [record(submission), record(preparation), submitted["script"]])
    checked.extend(record(directory / name) for name in (
        "plan.json", "execution_started.json", "execution.json", "run/complete.json",
        "run/started.json", "run/native/started.json", "run/score.json", "run/data.json"))
    validate(directory / "plan.json", PLAN_SHA)
    return dict(status="execution_verified_pending_scientific_admission", job_id=job,
                accounting=accounting, scheduler=scheduler, checked_records=checked,
                source=record(__file__), scientific_admission=False, controlled_timing=False,
                publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = verify(args.directory.resolve())
    save(args.output, result)
    print(result["status"])

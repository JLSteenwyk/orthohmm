"""One explicit full-review postprocessing attempt after successful diagnosis.

Execute the original reviewer in a fresh namespace; never inference or scoring.
The larger postprocessing allocation is not a change to native resource limits.
"""

import argparse
import json
import os
from pathlib import Path
import pwd
import re
import subprocess
import sys
import time

from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import ROOT, read, available_memory, validate_plan, validate_request
from benchmark_tools.run_review_gated_native_qfo_pairs import runtime_environment
from benchmark_tools.validate_native_factorial_outputs import require


REQUEST = ROOT / "benchmarks/work/native_factorial_launch_20261004/request_08_receipt_amended.json"
REQUEST_SHA = "35f8ef9b1d7abc1f574b8e7e3c55bc91c68c9f8f74247a9a816754a41a93fb6a"
REVIEWER = ROOT / "benchmark_tools/review_native_factorial_attempt.py"
REVIEWER_SHA = "63e7d7fdda52afa7a36eecd89492a2684e8310260febbfb5adf6d23f88fd26de"
VALIDATOR = ROOT / "benchmark_tools/validate_native_factorial_outputs.py"
VALIDATOR_SHA = "3357637503f35238f654edea5c4c12bd293f8d20908c141421a2d40278af7c8d"
DIAGNOSTIC = ROOT / "benchmarks/work/native_review_sigsegv_diagnostic_20261006_v1/outputs.json"
DIAGNOSTIC_SUBMISSION = ROOT / "benchmark_tools/results/native_review_sigsegv_submission_22734.json"
DIAGNOSTIC_SUBMISSION_SHA = "4f1df72adca42257ec48d8b3585a4e1f4f462d6aefb46098fd913df5c9db4d14"
PYTHON = ROOT / "benchmarks/work/native_factorial_review_py310_20261004/bin/python"
BATCH = ROOT / "benchmark_tools/results/fault_reported_native_review_22444_20261006.sh"
CONTROL = ROOT / "benchmarks/work/native_factorial_fault_reported_review_22444_v1"
DESTINATION = ROOT / "benchmarks/work/native_factorial_terminal_review_22444_fault_reported_v1"
MEMORY = 128 * 2**30


def diagnostic_scheduler_gate(scheduler):
    require(scheduler.get("JobIDRaw") == "22734" and scheduler.get("State") == "COMPLETED"
        and scheduler.get("ExitCode") == "0:0" and scheduler.get("NodeList") == "bizon"
        and scheduler.get("AllocCPUS") == "2"
        and scheduler.get("ReqMem") in {"32G", "32Gn", "32768M", "32768Mn"},
        "Require successful original diagnostic in its declared envelope")


def diagnostic_gate(output, scheduler, request_ref, plan_ref):
    diagnostic_scheduler_gate(scheduler)
    require(output.get("schema") == "native_factorial_output_review_v1"
        and output.get("status") == "native_outputs_validated" and output.get("job_id") == 22444
        and type(output.get("index")) is int and output["index"] == 8
        and output.get("cell") == "p0_c1_r0" and output.get("request") == request_ref
        and output.get("plan") == plan_ref and output.get("source") == record(VALIDATOR)
        and output.get("native_outputs_validated") is True
        and output.get("terminal_scheduler_confirmed") is True
        and all(output.get(k) is False for k in ("accuracy_evaluated", "terminal_reviewed",
            "resource_measurements_admitted", "next_identity_authorized", "uncontended_timing")),
        "Require bound original standalone semantic result, not scientific admission")


def allocation_gate(raw, job, digest):
    lines = [line for line in raw.splitlines() if line.strip()]
    require(len(lines) == 1, "Require one fresh full-review controller record")
    pairs = re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", lines[0])
    fields = dict(pairs)
    require(len(fields) == len(pairs), "Duplicate full-review controller field")
    expected = dict(JobId=str(job), JobState="RUNNING", Partition="gpu", NodeList="bizon",
        NumNodes="1", NumCPUs="2", NumTasks="1", MinMemoryNode="128G", TimeLimit="06:00:00",
        Requeue="0", Restarts="0", Command=str(BATCH), WorkDir=str(ROOT), Comment=digest,
        UserId=f"{pwd.getpwuid(os.getuid()).pw_name}({os.getuid()})")
    expected["CPUs/Task"] = "2"
    require(type(job) is int and job > 22734 and all(fields.get(k) == v for k, v in expected.items())
        and not {"ArrayJobId", "ArrayTaskId", "HetJobId", "HetJobOffset"}.intersection(fields),
        "Require a distinct owned two-CPU/128GiB full-review attempt")
    return fields


def execute(diagnostic_sha256, worker_sha256):
    source = record(__file__)
    require(source["sha256"] == worker_sha256, "Full-review wrapper source changed")
    request_ref = record(REQUEST)
    require(request_ref["sha256"] == REQUEST_SHA and record(REVIEWER)["sha256"] == REVIEWER_SHA
        and record(VALIDATOR)["sha256"] == VALIDATOR_SHA, "Original review/request sources changed")
    request = read(request_ref)
    plan = read(request["plan"])
    run = validate_plan(plan)[8]
    validate_request(request, request["plan"], 22444)
    require(run["cell"] == "p0_c1_r0", "Native cell differs")
    text, scheduler = accounting(22734, include_memory=True)
    diagnostic_scheduler_gate(scheduler)
    producer_ref = record(DIAGNOSTIC_SUBMISSION)
    require(producer_ref["sha256"] == DIAGNOSTIC_SUBMISSION_SHA, "Original diagnostic submission changed")
    producer = read(producer_ref)
    require(producer.get("schema") == "native_reviewer_signal11_diagnostic_submission_v1"
        and producer.get("job_id") == 22734 and producer.get("native_job_id") == 22444
        and producer.get("original_failed_reviewer") == 22445 and producer.get("index") == 8
        and producer.get("request") == request_ref and producer.get("destination") == str(DIAGNOSTIC.parent)
        and producer.get("automatic_retry") is False, "Diagnostic producer identity differs")
    check(producer["batch"])
    # Observe a future digest only after the pinned producer is terminal and successful.
    diagnostic_ref = record(DIAGNOSTIC)
    require(diagnostic_sha256 is None or diagnostic_ref["sha256"] == diagnostic_sha256,
        "Diagnostic result digest differs")
    diagnostic = read(diagnostic_ref)
    diagnostic_gate(diagnostic, scheduler, request_ref, request["plan"])
    failure_command = ["sacct", "-j", "22445", "--parsable2", "--noheader", "--format=JobIDRaw,State,ExitCode"]
    failed = subprocess.run(failure_command, capture_output=True, text=True, check=True)
    require("22445|FAILED|0:11" in failed.stdout.splitlines(), "Original reviewer failure not preserved")
    original = ROOT / "benchmarks/work/native_factorial_terminal_review_22444"
    require(not (original / "review.json").exists() and not (original / "failure.json").exists(),
        "Original failed-review namespace differs")
    refs = [source, record(BATCH), request_ref, request["plan"], diagnostic_ref, producer_ref, producer["batch"],
        record(REVIEWER), record(VALIDATOR),
        *diagnostic["checked_files"], *diagnostic["evidence"], *plan["helper_sources"], *plan["evidence"]]
    for ref in refs:
        check(ref)
    runtime = runtime_environment(dict(python_invocation_path=str(PYTHON), python_binary=record(PYTHON)))
    job = int(os.environ.get("SLURM_JOB_ID", "0"))
    require(os.environ.get("SLURM_CPUS_PER_TASK") == "2", "Require scheduled full review")
    query = ["scontrol", "show", "job", str(job), "--oneliner"]
    observed = subprocess.run(query, capture_output=True, text=True, check=True, timeout=5)
    allocation = allocation_gate(observed.stdout, job, worker_sha256 if diagnostic_sha256 is None else diagnostic_sha256)
    raw_meminfo = Path("/proc/meminfo").read_text()
    capacity = available_memory(raw_meminfo)
    require(capacity >= MEMORY, "Unsafe full-review available memory")
    for path in (CONTROL, DESTINATION):
        require(path.resolve() == path and not path.exists() and not path.is_symlink(),
            "Require fresh direct postprocessing namespaces")
    CONTROL.mkdir(parents=True, exist_ok=False)
    command = [str(PYTHON), "-B", "-X", "faulthandler", "-m", "benchmark_tools.review_native_factorial_attempt",
        "--request", str(REQUEST), "--request-sha256", REQUEST_SHA, "--output-directory", str(DESTINATION)]
    report = dict(schema="fault_reported_native_factorial_review_attempt_v1", status="running_original_full_review",
        source=source, diagnostic=diagnostic_ref, diagnostic_submission=producer_ref,
        diagnostic_digest_mode="observed_after_producer_completion" if diagnostic_sha256 is None else "supplied_after_producer_completion",
        diagnostic_accounting=text, diagnostic_scheduler=scheduler,
        job_id=job, native_job_id=22444, failed_reviewer_job_id=22445, request=request_ref, plan=request["plan"],
        reviewer=record(REVIEWER), command=command, destination=str(DESTINATION), runtime=runtime,
        scheduler=dict(command=query, stdout=observed.stdout, stderr=observed.stderr, fields=allocation),
        original_failure_accounting=dict(command=failure_command, stdout=failed.stdout, stderr=failed.stderr),
        raw_meminfo=raw_meminfo, available_memory_bytes=capacity, checked_records=refs,
        started_monotonic_ns=time.monotonic_ns(), original_full_review_reexecuted=True,
        native_inference_reexecuted=False, accuracy_evaluated=False, scientific_timings_admitted=False,
        next_identity_authorized=False, automatic_retry=False, publication_ready=False)
    save(CONTROL / "preflight.json", report)
    try:
        with (CONTROL / "stdout.txt").open("x") as stdout, (CONTROL / "stderr.txt").open("x") as stderr:
            child = subprocess.run(command, cwd=ROOT, stdout=stdout, stderr=stderr, check=False)
        report["child_exit_code"] = child.returncode
        require(child.returncode == 0, "Fault-reported original full reviewer failed; retain this attempt")
        review_ref = record(DESTINATION / "review.json")
        review = read(review_ref)
        require(review.get("schema") == "native_factorial_terminal_review_v1"
            and review.get("status") == "native_success" and review.get("job_id") == 22444
            and review.get("index") == 8 and review.get("request") == request_ref
            and review.get("source") == record(REVIEWER) and review.get("native_outputs_validated") is True
            and review.get("terminal_reviewed") is True and review.get("next_identity_authorized") is True,
            "Original full-review result differs")
        for ref in refs:
            check(ref)
        report.update(status="original_full_review_returned_pending_independent_completion_check", review=review_ref)
    except BaseException as error:
        report.update(status="fault_reported_full_review_failed_retained", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_monotonic_ns"] = time.monotonic_ns()
        report["limitations"] = [
            "Explicit new postprocessing attempt; original reviewer/downstream failures remain unchanged.",
            "128GiB is a larger review-only limit, not the native scientific resource envelope or a proven crash fix.",
            "The original full reviewer repeats every gate; standalone diagnosis is not substituted for admission.",
            "Producer completion and independent review remain required before conversion or next-identity release.",
            "Shared-host postprocessing timings have unknown, potentially tool-dependent contention effects.",
            "No inference, scoring, automatic retry, source patching, monkeypatching or isolated-performance claim."]
        save(CONTROL / "results.json", report)
    return record(CONTROL / "results.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostic-sha256",
        help="Optional known digest; otherwise observe it only after successful producer completion")
    parser.add_argument("--worker-sha256", required=True)
    args = parser.parse_args()
    print(json.dumps(execute(args.diagnostic_sha256, args.worker_sha256), sort_keys=True))

"""One prospectively gated full-review attempt after native11 standalone diagnosis."""

import argparse
import json
import os
from pathlib import Path
import pwd
import re
import subprocess
import time

from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting
from benchmark_tools.native_factorial_allocated_execution import amendment, validate_request
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import ROOT, available_memory, read
from benchmark_tools.run_review_gated_native_qfo_pairs import runtime_environment
from benchmark_tools.validate_native_factorial_outputs import require


REQUEST = ROOT / "benchmarks/work/native_factorial_launch_20261004/request_11_allocated_v1.json"
REQUEST_SHA = "7bf63b80bd5932b9edbd1b2c5ff3fb77f5557e4f6c64e077045d6a50c8d366a1"
REVIEWER = ROOT / "benchmark_tools/review_allocated_native_factorial_attempt.py"
REVIEWER_SHA = "c8707fdb1855a0d73fba537c001b4a355faff094a38a177a05e179148924c929"
VALIDATOR = ROOT / "benchmark_tools/validate_allocated_native_factorial_outputs.py"
VALIDATOR_SHA = "146f1333715b5ef52ee17de855fad584d1efc91c7a1154edcd269829b5a871be"
SEMANTIC = ROOT / "benchmark_tools/validate_native_factorial_outputs.py"
SEMANTIC_SHA = "3357637503f35238f654edea5c4c12bd293f8d20908c141421a2d40278af7c8d"
PREPARATION = ROOT / "benchmark_tools/results/native11_standalone_diagnostic_preparation_20261008_v1.json"
PREPARATION_SHA = "fc2ffe1c6b46647648fa4c1e8068e260bc605bb7328f0f2e1946090162174f98"
SUBMISSION = ROOT / "benchmark_tools/results/native11_standalone_diagnostic_held_20261008_v1.json"
SUBMISSION_SHA = "95da0070b2e475e4c0885fdcc25ffb48f4cced939302bcf67731d1b116958e77"
DIAGNOSTIC = ROOT / "benchmarks/work/native11_standalone_diagnostic_20261008_v1/outputs.json"
PYTHON = ROOT / "benchmarks/work/native_factorial_review_py310_20261004/bin/python"
BATCH = ROOT / "benchmark_tools/results/native11_fault_reported_review_20261008_v1.sh"
CONTROL = ROOT / "benchmarks/work/native11_fault_reported_review_20261008_v1"
DESTINATION = ROOT / "benchmarks/work/native11_full_terminal_review_20261008_v1"


def diagnostic_scheduler_gate(fields):
    require(fields.get("JobIDRaw") == "24031" and fields.get("State") == "COMPLETED"
        and fields.get("ExitCode") == "0:0" and fields.get("NodeList") == "bizon"
        and fields.get("AllocCPUS") == "2"
        and fields.get("ReqMem") in {"32G", "32Gn", "32768M", "32768Mn"},
        "Require successful diagnostic24031 in its declared resource envelope")


def diagnostic_gate(output, fields, request_ref, plan_ref, amendment_ref, validator_ref):
    diagnostic_scheduler_gate(fields)
    require(output.get("schema") == "allocated_native_factorial_output_review_v1"
        and output.get("status") == "native_outputs_validated" and output.get("job_id") == 23985
        and type(output.get("index")) is int and output["index"] == 11
        and output.get("cell") == "p1_c1_r0" and output.get("request") == request_ref
        and output.get("plan") == plan_ref and output.get("amendment") == amendment_ref
        and output.get("source") == validator_ref and output.get("native_outputs_validated") is True
        and output.get("terminal_scheduler_confirmed") is True
        and all(output.get(k) is False for k in ("accuracy_evaluated", "terminal_reviewed",
            "resource_measurements_admitted", "next_identity_authorized", "uncontended_timing")),
        "Require bound standalone semantic result, not full-review admission")


def allocation_gate(raw, job, digest):
    lines = [line for line in raw.splitlines() if line.strip()]
    require(len(lines) == 1, "Require one fresh full-review controller record")
    pairs = re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", lines[0])
    fields = dict(pairs)
    require(len(fields) == len(pairs), "Duplicate full-review controller field")
    expected = dict(JobId=str(job), JobName="ohmm_native11_fullreview", JobState="RUNNING",
        Partition="gpu", NodeList="bizon", NumNodes="1", NumCPUs="2", NumTasks="1",
        MinMemoryNode="128G", TimeLimit="06:00:00", Requeue="0", Restarts="0",
        Command=str(BATCH), WorkDir=str(ROOT), Comment=digest,
        UserId=f"{pwd.getpwuid(os.getuid()).pw_name}({os.getuid()})")
    expected["CPUs/Task"] = "2"
    require(type(job) is int and job > 24031
        and all(fields.get(k) == v for k, v in expected.items())
        and not {"ArrayJobId", "ArrayTaskId", "HetJobId", "HetJobOffset"}.intersection(fields),
        "Require a distinct owned two-CPU/128GiB full-review job")
    return fields


def execute(diagnostic_sha256, worker_sha256):
    source = record(__file__)
    require(source["sha256"] == worker_sha256, "Full-review wrapper changed")
    request_ref = record(REQUEST)
    require(request_ref["sha256"] == REQUEST_SHA, "Original request changed")
    for path, digest in ((REVIEWER, REVIEWER_SHA), (VALIDATOR, VALIDATOR_SHA), (SEMANTIC, SEMANTIC_SHA),
                         (PREPARATION, PREPARATION_SHA), (SUBMISSION, SUBMISSION_SHA)):
        require(record(path)["sha256"] == digest, "Original source or diagnostic binding changed")
    request = read(request_ref)
    execution, plan = amendment(request["amendment"])
    validate_request(request, request["amendment"], execution, 23985)
    require(plan["runs"][11]["cell"] == "p1_c1_r0", "Native cell differs")
    text, fields = accounting(24031, include_memory=True)
    diagnostic_scheduler_gate(fields)
    submission_ref = record(SUBMISSION)
    submission = read(submission_ref)
    require(submission.get("job_id") == 24031 and submission.get("held_comparison_passed") is True
        and submission.get("submission_count") == 1 and submission.get("automatic_retry") is False,
        "Diagnostic submission differs")
    check(submission["batch"])
    # Only observe the result after its original producer is successfully terminal.
    diagnostic_ref = record(DIAGNOSTIC)
    require(diagnostic_ref["sha256"] == diagnostic_sha256, "Diagnostic result checksum differs")
    diagnostic = read(diagnostic_ref)
    diagnostic_gate(diagnostic, fields, request_ref, request["plan"], request["amendment"], record(VALIDATOR))
    preparation_ref = record(PREPARATION)
    preparation = read(preparation_ref)
    failure_command = ["sacct", "-X", "-j", "23986", "-n", "-P", "--format=JobIDRaw,State,ExitCode"]
    failed = subprocess.run(failure_command, capture_output=True, text=True, check=True, timeout=10)
    require(failed.stdout.strip() == "23986|FAILED|0:11", "Original review outcome differs")
    original = ROOT / "benchmarks/work/allocated_native_factorial_terminal_review_23985_v1"
    require({p.name for p in original.iterdir()} ==
        {Path(ref["path"]).name for ref in preparation["original_partial_review_files"]},
        "Original failed-review inventory differs")
    refs = [source, record(BATCH), request_ref, request["plan"], request["amendment"],
        record(REVIEWER), record(VALIDATOR), record(SEMANTIC), preparation_ref, submission_ref,
        submission["batch"], diagnostic_ref, preparation["original_review_log"],
        *preparation["original_partial_review_files"], *diagnostic["evidence"], *diagnostic["checked_files"]]
    for ref in refs:
        check(ref)
    runtime = runtime_environment(dict(python_invocation_path=str(PYTHON), python_binary=record(PYTHON)))
    job = int(os.environ.get("SLURM_JOB_ID", "0"))
    require(os.environ.get("SLURM_CPUS_PER_TASK") == "2", "Require scheduled two-CPU full review")
    query = ["scontrol", "show", "job", str(job), "--oneliner"]
    observed = subprocess.run(query, capture_output=True, text=True, check=True, timeout=10)
    allocation = allocation_gate(observed.stdout, job, diagnostic_sha256)
    raw_memory = Path("/proc/meminfo").read_text()
    capacity = available_memory(raw_memory)
    require(capacity >= 128*2**30, "Unsafe review capacity")
    for path in (CONTROL, DESTINATION):
        require(path.resolve() == path and not path.exists() and not path.is_symlink(), "Full-review namespace used")
    CONTROL.mkdir(exist_ok=False)
    command = [str(PYTHON), "-B", "-X", "faulthandler", "-m", "benchmark_tools.review_allocated_native_factorial_attempt",
        "--request", str(REQUEST), "--request-sha256", REQUEST_SHA, "--output-directory", str(DESTINATION)]
    report = dict(schema="native11_fault_reported_full_review_attempt_v1", status="running_original_full_review",
        source=source, batch=record(BATCH), request=request_ref, diagnostic=diagnostic_ref,
        diagnostic_submission=submission_ref, diagnostic_accounting=text, diagnostic_scheduler=fields,
        job_id=job, native_job_id=23985, failed_reviewer_job_id=23986, command=command,
        runtime=runtime, checked_records=refs, destination=str(DESTINATION),
        scheduler=dict(command=query, stdout=observed.stdout, stderr=observed.stderr, fields=allocation),
        original_failure_accounting=dict(command=failure_command, stdout=failed.stdout, stderr=failed.stderr),
        available_memory_bytes=capacity, raw_meminfo=raw_memory, started_monotonic_ns=time.monotonic_ns(),
        original_full_review_reexecuted=True, native_inference_reexecuted=False, automatic_retry=False,
        historical_failure_cause_established=False, accuracy_evaluated=False,
        scientific_timings_admitted=False, next_identity_authorized=False, publication_ready=False)
    save(CONTROL / "preflight.json", report)
    try:
        with (CONTROL / "stdout.txt").open("x") as stdout, (CONTROL / "stderr.txt").open("x") as stderr:
            child = subprocess.run(command, cwd=ROOT, stdout=stdout, stderr=stderr, check=False)
        report["child_exit_code"] = child.returncode
        require(child.returncode == 0, "Original full reviewer failed; preserve this attempt without retry")
        review_ref = record(DESTINATION / "review.json")
        review = read(review_ref)
        require(review.get("schema") == "allocated_native_factorial_terminal_review_v1"
            and review.get("status") == "native_success" and review.get("job_id") == 23985
            and review.get("index") == 11 and review.get("request") == request_ref
            and review.get("source") == record(REVIEWER)
            and all(review.get(k) is True for k in ("native_outputs_validated", "terminal_reviewed",
                "next_identity_authorized", "primary_resources_replayed", "shared_host_resources_reviewed")),
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
            "One explicit prospective postprocessing attempt; original23986 failure remains retained.",
            "The original full reviewer repeats every gate; standalone output is not substituted for admission.",
            "Producer completion and independent review are still required for dependent conversion or native12.",
            "128GiB review envelope is unchanged from23986 and is not proof of an OOM or crash fix.",
            "Shared-host postprocessing timing has unknown, potentially tool-dependent contention effects.",
            "No inference, scoring, GC alteration, source patch, monkeypatch or automatic retry."]
        save(CONTROL / "results.json", report)
    return record(CONTROL / "results.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostic-sha256", required=True)
    parser.add_argument("--worker-sha256", required=True)
    args = parser.parse_args()
    print(json.dumps(execute(args.diagnostic_sha256, args.worker_sha256), sort_keys=True))

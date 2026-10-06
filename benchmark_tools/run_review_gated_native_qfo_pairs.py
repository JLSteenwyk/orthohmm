"""Run unchanged native QfO conversion only after its bound reviewer completes."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time

from benchmark_tools import prepare_native_factorial_qfo_pairs as conversion
from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import ROOT, available_memory, read
from benchmark_tools.validate_native_factorial_outputs import require


CONVERTER_SHA = "472f26a140cf89aed36f8212a31fef7d8f5061204dac2e0ce80c3c7e6746630e"


def submission_binding(submission_ref):
    submission = read(submission_ref)
    require(submission.get("schema") == "native_factorial_review_held_submission_v1"
        and type(submission.get("job_id")) is int and submission["job_id"] > 0
        and type(submission.get("native_job_id")) is int and submission["native_job_id"] > 0
        and submission["job_id"] != submission["native_job_id"]
        and type(submission.get("index")) is int and 6 <= submission["index"] < 13
        and submission.get("automatic_retry") is False
        and submission.get("accuracy_evaluated") is False,
        "Require an explicit bound native-QfO reviewer submission")
    request_ref = submission["request"]
    request = read(request_ref)
    conversion.validate_request(request, request["plan"], submission["native_job_id"])
    run = conversion.validate_plan(read(request["plan"]))[submission["index"]]
    require(request["index"] == submission["index"] == run["index"], "Reviewer/native index differs")
    conversion.conversion_kind(run)
    expected = ROOT / "benchmarks/work" / f'native_factorial_terminal_review_{submission["native_job_id"]}'
    require(submission.get("destination") == str(expected), "Unexpected reviewer output location")
    require(submission["review_source"] == record(ROOT / "benchmark_tools/review_native_factorial_attempt.py"),
            "Bound reviewer source changed")
    check(submission["batch"])
    return submission, request_ref, request, run


def reviewer_completed(fields, submission):
    require(fields.get("JobIDRaw") == str(submission["job_id"])
        and fields.get("State") == "COMPLETED" and fields.get("ExitCode") == "0:0"
        and fields.get("NodeList") == "bizon" and fields.get("AllocCPUS") == "2"
        and fields.get("ReqMem") in {"32G", "32Gn", "32768M", "32768Mn"},
        "Reviewer is not successfully completed in its bound resource envelope")


def runtime_environment(submission):
    invocation = submission["python_invocation_path"]
    require(sys.executable == invocation and sys.prefix == str(Path(invocation).parent.parent)
        and record(invocation) == submission["python_binary"]
        and sys.version_info[:3] == (3, 10, 13), "Require the original Python3.10 venv invocation")
    import Bio
    import numpy
    import psutil
    versions = dict(Bio=Bio.__version__, numpy=numpy.__version__, psutil=psutil.__version__)
    require(versions == dict(Bio="1.87", numpy="2.2.6", psutil="7.2.2"), "Review environment versions changed")
    return dict(invocation=invocation, prefix=sys.prefix, version=sys.version,
                binary=record(invocation), packages=versions)


def capacity():
    meminfo = Path("/proc/meminfo").read_text()
    return dict(observed_unix_ns=time.time_ns(), available_memory_bytes=available_memory(meminfo),
                available_disk_bytes=shutil.disk_usage(ROOT / "benchmarks/work").free)


def execute(submission_ref, worker_sha256):
    source = record(__file__)
    require(source["sha256"] == worker_sha256, "Gated-conversion worker source changed")
    submission, request_ref, request, run = submission_binding(submission_ref)
    native_job = submission["native_job_id"]
    parent = ROOT / "benchmarks/work"
    directory = parent / f"native_factorial_qfo_conversion_gate_{native_job}"
    destination = parent / f"native_factorial_qfo_pairs_{native_job}"
    for path in (directory, destination):
        require(path.is_absolute() and path.resolve() == path, "Require direct output namespaces")
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    require(os.environ.get("SLURM_CPUS_PER_TASK") == "2"
        and os.environ.get("SLURM_JOB_ID", "").isdigit()
        and int(os.environ["SLURM_JOB_ID"]) not in {native_job, submission["job_id"]},
        "Require a separate scheduled two-CPU conversion job")
    directory.mkdir(parents=True, exist_ok=False)
    report = dict(schema="review_gated_native_qfo_conversion_v1", status="validating_after_reviewer",
        source=source, reviewer_submission=submission_ref, request=request_ref,
        job_id=os.environ["SLURM_JOB_ID"], native_job_id=native_job, reviewer_job_id=submission["job_id"],
        index=run["index"], cell=run["cell"], destination=str(destination),
        started_monotonic_ns=time.monotonic_ns(), native_inference_reexecuted=False,
        accuracy_evaluated=False, next_identity_authorized=False, automatic_retry=False,
        publication_ready=False, conversion_started=False)
    try:
        text, fields = accounting(submission["job_id"], include_memory=True)
        reviewer_completed(fields, submission)
        report.update(reviewer_accounting=text, reviewer_scheduler=fields)
        report["runtime"] = runtime_environment(submission)
        review_ref = record(Path(submission["destination"]) / "review.json")
        review = read(review_ref)
        kind = conversion.admit_conversion(review, request_ref, request, run)
        converter_ref = record(conversion.__file__)
        require(converter_ref["sha256"] == CONVERTER_SHA, "Original conversion source changed")
        context = capacity()
        report.update(terminal_review=review_ref, converter=converter_ref, conversion_kind=kind,
                      resource_context=context)
        require(context["available_memory_bytes"] >= 32 * 2**30
            and context["available_disk_bytes"] >= 128 * 2**30, "Unsafe conversion memory/disk capacity")
        for ref in (source, submission_ref, request_ref, review_ref, converter_ref,
                    submission["review_source"], submission["batch"]):
            check(ref)
        report.update(status="running_unchanged_converter", conversion_started=True)
        save(directory / "preflight.json", report)
        converted = conversion.prepare(request_ref, review_ref, destination)
        stage = read(converted)
        require(converted["path"] == str(destination / "results.json")
            and stage.get("status") == "full_native_factorial_qfo_pairs_prepared_unscored"
            and stage.get("native_job_id") == native_job and stage.get("native_index") == run["index"]
            and stage.get("cell") == run["cell"] and stage.get("job_id") == report["job_id"]
            and stage.get("source") == converter_ref and stage.get("terminal_review") == review_ref
            and stage.get("accuracy_evaluated") is False and stage.get("next_identity_authorized") is False,
            "Unexpected converted stage identity or scope")
        for ref in (source, submission_ref, request_ref, review_ref, converter_ref):
            check(ref)
        report.update(status="review_gated_native_qfo_pairs_prepared_unscored", conversion=converted)
    except BaseException as error:
        report.update(status="review_gated_native_qfo_conversion_failed_retained",
                      error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_monotonic_ns"] = time.monotonic_ns()
        report["limitations"] = [
            "Postprocessing only: endpoint assessment and independent scientific admission remain required.",
            "A completed reviewer alone is insufficient; the unchanged native conversion admission gates must pass.",
            "Future review digest is observed only after producer completion and validated against the pinned request/source/output identities.",
            "Shared-host capacity observation does not establish isolation or budget adequacy for the whole conversion.",
            "No automatic inference, recovery, scoring, retry, successor release or timing repair."]
        save(directory / "results.json", report)
    return record(directory / "results.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review-submission", required=True, type=Path)
    parser.add_argument("--review-submission-sha256", required=True)
    parser.add_argument("--worker-sha256", required=True)
    args = parser.parse_args()
    submission_ref = record(args.review_submission)
    require(submission_ref["sha256"] == args.review_submission_sha256, "Reviewer submission checksum differs")
    print(json.dumps(execute(submission_ref, args.worker_sha256), sort_keys=True))

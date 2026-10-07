"""One prospective terminal review after pinned coverage/context diagnosis; no inference retry."""

import argparse
import json
import os
from pathlib import Path
import time
import traceback

from benchmark_tools.diagnose_native_partition import Bindings


ROOT = Path(__file__).resolve().parents[1]
CONTROL = ROOT / "benchmarks/work/native10_review_finalization_20261007_v1"
PARTITION = ROOT / "benchmark_tools/results/native10_partition_diagnosis_20261007_v1.json"
PARTITION_SHA = "4d8ce7b56b07d9a14cc5ee511ca13ac7021cc865908b28fcd5db33c85f4bc4e1"
CONTEXT = ROOT / "benchmark_tools/results/native10_semantic_context_diagnosis_20261007_v1.json"
CONTEXT_SHA = "d54f2062f1dd37a25a69c03b6e4993466f96ac2c5313692bc60f4580a732ab1c"
FAILURE_SHA = "bad4de5092098cf461ad618a3b9527cf2f02800d35a09e37fe0e413447560e99"
REQUEST_SHA = "1355ae3b73a9c8496133bd25a6cbf136d149598ccab680672e81f86482338910"


def decision(partition, context, failure, request, job, cpus):
    if (type(job) is not int or job <= 0 or job in (23902, 23910) or cpus != "2"
            or partition.get("index") != 10 or partition.get("job_id") != 23902
            or partition.get("status") != "diagnosis_completed"
            or partition.get("frozen_root_coverage_gate", {}).get("status") != "passed"
            or partition.get("all_output_partitions_identical") is not True
            or partition.get("source_payload_matches_native_input") is not True
            or partition.get("parsers_identical") is not True
            or not partition.get("stages")
            or any(s.get("complete_unique_input_partition") is not True for s in partition["stages"].values())
            or context.get("index") != 10 or context.get("job_id") != 23902
            or context.get("status") != "semantic_probe_passed"
            or context.get("semantic_result", {}).get("native_outputs_validated") is not True
            or failure.get("index") != 10 or failure.get("job_id") != 23902
            or failure.get("status") != "terminal_factorial_review_failed"
            or request.get("index") != 10 or request.get("job_id") != 23902):
        raise ValueError("Require both successful pinned diagnoses and a distinct scheduled two-CPU review")
    for value in (partition, context, failure):
        if any(value.get(k) is not False for k in (
                "accuracy_evaluated", "terminal_reviewed", "next_identity_authorized", "automatic_retry")):
            raise ValueError("Diagnosis or original failure was promoted to admission")
    return dict(status="one_shot_fresh_review_permitted", native_job_id=23902,
                failed_review_job_id=23910, new_review_job_id=job, native_retry=False,
                original_failure_retained=True, automatic_retry=False,
                historical_failure_cause_established=False)


def finalize(job, cpus):
    from benchmark_tools.review_allocated_native_factorial_attempt import review
    from benchmark_tools.run_native_factorial_cost import available_memory

    if CONTROL.exists() or CONTROL.is_symlink():
        raise ValueError("Finalization identity already used; no retry or overwrite")
    evidence = Bindings()
    partition = evidence.read(PARTITION)
    context = evidence.read(CONTEXT)
    if (evidence.files[str(PARTITION)]["sha256"] != PARTITION_SHA
            or evidence.files[str(CONTEXT)]["sha256"] != CONTEXT_SHA):
        raise ValueError("Pinned diagnostic checksum differs")
    for ref in [*partition["evidence"], *context["evidence"]]:
        evidence.bind(ref["path"], ref)
    failure_ref = partition["original_failure"]
    failure = evidence.read(failure_ref["path"], failure_ref)
    evidence.bind(failure["source"]["path"], failure["source"])
    request_ref = failure["request"]
    request = evidence.read(request_ref["path"], request_ref)
    if failure_ref["sha256"] != FAILURE_SHA or request_ref["sha256"] != REQUEST_SHA:
        raise ValueError("Original failure/request checksum differs")
    permitted = decision(partition, context, failure, request, job, cpus)
    raw_memory = Path("/proc/meminfo").read_text()
    capacity = available_memory(raw_memory)
    if capacity < 128 * 2**30:
        raise ValueError("Unsafe review capacity; no automatic retry")
    evidence.bind(Path(__file__).resolve())
    checked = evidence.finish()
    CONTROL.mkdir(exist_ok=False)
    preflight = dict(permitted, observed_unix_ns=time.time_ns(), evidence=checked,
                     available_memory_bytes=capacity, minimum_available_memory_bytes=128 * 2**30,
                     raw_meminfo=raw_memory, scientific_outputs_admitted=False)
    (CONTROL / "preflight.json").write_text(json.dumps(preflight, indent=2, sort_keys=True) + "\n")
    try:
        # The unchanged original reviewer performs every check and writes its own source identity.
        result = review(request_ref, CONTROL / "review")
    except BaseException as error:
        (CONTROL / "finalization_failure.json").write_text(json.dumps(dict(
            status="fresh_review_failed", new_review_job_id=job, native_job_id=23902,
            original_failure=failure_ref, automatic_retry=False, accuracy_evaluated=False,
            error_type=type(error).__name__, error=str(error),
            traceback="".join(traceback.format_exception(type(error), error, error.__traceback__))),
            indent=2, sort_keys=True) + "\n")
        raise
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    job = os.environ.get("SLURM_JOB_ID", "")
    if not job.isdecimal():
        raise ValueError("Require a scheduled review")
    print(json.dumps(finalize(int(job), os.environ.get("SLURM_CPUS_PER_TASK")), sort_keys=True))

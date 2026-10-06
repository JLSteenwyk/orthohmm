"""Run unchanged native QfO assessment after its bound conversion completes."""

import argparse
import json
import os
from pathlib import Path
import time

from benchmark_tools import run_native_factorial_qfo_assessment as assessment
from benchmark_tools import run_review_gated_native_qfo_pairs as conversion_gate
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import ROOT, read
from benchmark_tools.validate_native_factorial_outputs import require


ASSESSMENT_SHA = "ba051efa4b2dd9f72a4464f56a891b4742fa3a1c6999f8b95a65d0805c77dafe"


def submission_binding(submission_ref):
    submission = read(submission_ref)
    require(submission.get("schema") == "review_gated_native_qfo_conversion_held_submission_v1"
        and type(submission.get("job_id")) is int and submission["job_id"] > 0
        and type(submission.get("native_job_id")) is int and submission["native_job_id"] > 0
        and type(submission.get("reviewer_job_id")) is int and submission["reviewer_job_id"] > 0
        and type(submission.get("index")) is int and 6 <= submission["index"] < 13
        and submission.get("held_inspection_passed") is True
        and all(submission.get(k) is False for k in ("pair_conversion_completed", "accuracy_evaluated",
            "automatic_retry", "next_identity_authorized", "publication_ready")),
        "Require a bound inspected conversion submission, not an assumed completed stage")
    reviewer, request_ref, request, run = conversion_gate.submission_binding(submission["reviewer_submission"])
    require(submission["native_job_id"] == reviewer["native_job_id"]
        and submission["reviewer_job_id"] == reviewer["job_id"]
        and submission["index"] == run["index"] and submission["cell"] == run["cell"]
        and submission["request"] == request_ref and submission["plan"] == request["plan"]
        and submission["job_id"] not in {reviewer["job_id"], reviewer["native_job_id"]},
        "Conversion submission/native reviewer binding differs")
    require(submission["worker"] == record(conversion_gate.__file__)
        and submission["converter"] == record(conversion_gate.conversion.__file__)
        and submission["converter"]["sha256"] == conversion_gate.CONVERTER_SHA,
        "Conversion worker or original converter changed")
    check(submission["batch"])
    expected = [ROOT / "benchmarks/work" / f'{prefix}_{submission["native_job_id"]}'
        for prefix in ("native_factorial_qfo_conversion_gate", "native_factorial_qfo_pairs")]
    require(submission["output_namespaces"] == [str(p) for p in expected], "Unexpected conversion namespaces")
    return submission, reviewer, request, run, expected


def completed_conversion(fields, submission):
    require(fields.get("JobIDRaw") == str(submission["job_id"])
        and fields.get("State") == "COMPLETED" and fields.get("ExitCode") == "0:0"
        and fields.get("NodeList") == "bizon" and fields.get("AllocCPUS") == "2"
        and fields.get("ReqMem") in {"32G", "32Gn", "32768M", "32768Mn"},
        "Conversion is not successfully completed in its bound resource envelope")


def completed_gate(gate, submission, pairs_ref):
    require(gate.get("schema") == "review_gated_native_qfo_conversion_v1"
        and gate.get("status") == "review_gated_native_qfo_pairs_prepared_unscored"
        and gate.get("source") == submission["worker"] and gate.get("converter") == submission["converter"]
        and gate.get("reviewer_submission") == submission["reviewer_submission"]
        and gate.get("request") == submission["request"] and gate.get("job_id") == str(submission["job_id"])
        and gate.get("native_job_id") == submission["native_job_id"]
        and gate.get("reviewer_job_id") == submission["reviewer_job_id"]
        and gate.get("index") == submission["index"] and gate.get("cell") == submission["cell"]
        and gate.get("destination") == submission["output_namespaces"][1]
        and gate.get("conversion") == pairs_ref and gate.get("conversion_started") is True
        and all(gate.get(k) is False for k in ("native_inference_reexecuted", "accuracy_evaluated",
            "next_identity_authorized", "automatic_retry", "publication_ready")),
        "Require the successful bound conversion gate, not merely converter output")


def output_paths(run, native_job):
    return dict(gate=ROOT / "benchmarks/work" / f"native_factorial_qfo_assessment_gate_{native_job}",
        cwd=ROOT / "benchmarks/results/full_native_qfo_assessment_v1" / run["cell"],
        work=ROOT / "qfo_benchmark/w" / f'nq{run["index"]:02d}',
        results=ROOT / "qfo_benchmark/scoring" / f'full_native_{run["index"]:02d}')


def execute(submission_ref, worker_sha256):
    source = record(__file__)
    require(source["sha256"] == worker_sha256, "Assessment gate worker changed")
    submission, reviewer, request, run, conversion_paths = submission_binding(submission_ref)
    paths = output_paths(run, submission["native_job_id"])
    for path in paths.values():
        require(path.is_absolute() and path.resolve() == path, "Require direct assessment namespaces")
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    require(os.environ.get("SLURM_CPUS_PER_TASK") == "8"
        and os.environ.get("SLURM_JOB_ID", "").isdigit()
        and int(os.environ["SLURM_JOB_ID"]) not in {submission["job_id"], submission["native_job_id"],
                                                   submission["reviewer_job_id"]},
        "Require a separate scheduled eight-CPU assessment")
    directory = paths["gate"]
    directory.mkdir(parents=True, exist_ok=False)
    report = dict(schema="conversion_gated_native_qfo_assessment_v1", status="validating_after_conversion",
        source=source, conversion_submission=submission_ref, job_id=os.environ["SLURM_JOB_ID"],
        native_job_id=submission["native_job_id"], conversion_job_id=submission["job_id"],
        reviewer_job_id=submission["reviewer_job_id"], index=run["index"], cell=run["cell"],
        output_namespaces={k: str(p) for k, p in paths.items()}, started_monotonic_ns=time.monotonic_ns(),
        assessment_driver_invoked=False, endpoint_process_completed=False, accuracy_admitted=False,
        native_inference_reexecuted=False, automatic_retry=False, next_identity_authorized=False,
        publication_ready=False)
    try:
        text, scheduler = assessment.accounting(submission["job_id"], include_memory=True)
        completed_conversion(scheduler, submission)
        report.update(conversion_accounting=text, conversion_scheduler=scheduler,
                      runtime=conversion_gate.runtime_environment(reviewer))
        gate_ref, pairs_ref = [record(p / "results.json") for p in conversion_paths]
        gate, stage = read(gate_ref), read(pairs_ref)
        completed_gate(gate, submission, pairs_ref)
        assessment.validate_stage(stage, scheduler)
        require(stage["request"] == submission["request"] and stage["plan"] == submission["plan"]
            and stage["native_job_id"] == submission["native_job_id"]
            and stage["native_index"] == run["index"] and stage["cell"] == run["cell"]
            and stage["source"] == submission["converter"] and stage["terminal_review"] == gate["terminal_review"],
            "Converted stage differs from successful gate/request")
        driver_ref = record(assessment.__file__)
        require(driver_ref["sha256"] == ASSESSMENT_SHA, "Original assessment driver changed")
        context = conversion_gate.capacity()
        require(context["available_memory_bytes"] >= 64 * 2**30
            and context["available_disk_bytes"] >= 128 * 2**30, "Unsafe assessment memory/disk capacity")
        refs = [source, submission_ref, gate_ref, pairs_ref, driver_ref, submission["worker"],
                submission["batch"], submission["reviewer_submission"], submission["request"], submission["plan"]]
        for ref in refs:
            check(ref)
        report.update(status="invoking_unchanged_assessment", conversion_gate=gate_ref, pairs=pairs_ref,
                      driver=driver_ref, resource_context=context, assessment_driver_invoked=True)
        save(directory / "preflight.json", report)
        execution = assessment.run(ROOT, pairs_ref, str(submission["job_id"]))
        execution_ref = record(paths["cwd"] / "results.json")
        require(read(execution_ref) == execution
            and execution.get("status") == "process_succeeded_pending_independent_admission"
            and execution.get("job_id") == report["job_id"] and execution.get("source") == driver_ref
            and execution.get("pairs_manifest") == pairs_ref and execution.get("exit_code") == 0
            and execution.get("native_index") == run["index"] and execution.get("cell") == run["cell"]
            and execution.get("native_job_id") == submission["native_job_id"]
            and all(execution.get(k) == str(paths[k]) for k in ("cwd", "work", "results"))
            and all(execution.get(k) is False for k in ("accuracy_admitted", "publication_ready",
                "automatic_retry", "native_inference_reexecuted", "next_identity_authorized")),
            "Unexpected original assessment result identity or scope")
        for ref in refs:
            check(ref)
        report.update(status="native_qfo_assessment_process_succeeded_pending_independent_admission",
                      execution=execution_ref, endpoint_process_completed=True)
    except BaseException as error:
        report.update(status="conversion_gated_native_qfo_assessment_failed_retained",
                      error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_monotonic_ns"] = time.monotonic_ns()
        report["limitations"] = [
            "Original assessment rechecks native terminal review, outputs, mapping, all helper pins and scorer environment before endpoint execution.",
            "Successful endpoint execution is not independent scientific admission or publication readiness.",
            "Future conversion digests are observed only after producer completion and bound to exact submitted request/gate identities.",
            "Shared-host capacity observation does not establish isolation or whole-job budget adequacy.",
            "No inference, automatic retry, successor release, scientific admission or timing repair."]
        save(directory / "results.json", report)
    return record(directory / "results.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conversion-submission", type=Path, required=True)
    parser.add_argument("--conversion-submission-sha256", required=True)
    parser.add_argument("--worker-sha256", required=True)
    args = parser.parse_args()
    ref = record(args.conversion_submission)
    require(ref["sha256"] == args.conversion_submission_sha256, "Conversion submission checksum differs")
    print(json.dumps(execute(ref, args.worker_sha256), sort_keys=True))

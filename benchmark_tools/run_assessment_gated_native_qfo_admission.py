"""Run unchanged independent QfO admission after bound assessment completion."""

import argparse
import json
import os
from pathlib import Path
import time

from benchmark_tools import admit_native_factorial_qfo_assessment as admission
from benchmark_tools import run_conversion_gated_native_qfo_assessment as assessment_gate
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import ROOT, read
from benchmark_tools.validate_native_factorial_outputs import require


ADMISSION_SHA = "1e13da48cb3879e1ce46c041b9e2579eb69c845a8208560e7a959a82dab0f8eb"


def submission_binding(submission_ref):
    sub = read(submission_ref)
    require(sub.get("schema") == "conversion_gated_native_qfo_assessment_held_submission_v1"
        and all(type(sub.get(k)) is int and sub[k] > 0 for k in
            ("job_id", "native_job_id", "reviewer_job_id", "conversion_job_id"))
        and type(sub.get("index")) is int and 6 <= sub["index"] < 13
        and len({sub[k] for k in ("job_id", "native_job_id", "reviewer_job_id", "conversion_job_id")}) == 4
        and sub.get("held_inspection_passed") is True
        and all(sub.get(k) is False for k in ("assessment_completed", "accuracy_admitted",
            "automatic_retry", "next_identity_authorized", "publication_ready")),
        "Require a bound inspected assessment submission")
    conversion, reviewer, request, run, conversion_paths = assessment_gate.submission_binding(sub["conversion_submission"])
    paths = assessment_gate.output_paths(run, sub["native_job_id"])
    require(sub["native_job_id"] == conversion["native_job_id"]
        and sub["reviewer_job_id"] == conversion["reviewer_job_id"]
        and sub["conversion_job_id"] == conversion["job_id"] and sub["index"] == run["index"]
        and sub["cell"] == run["cell"] and sub["request"] == conversion["request"]
        and sub["plan"] == conversion["plan"] and sub["reviewer_submission"] == conversion["reviewer_submission"]
        and sub["output_namespaces"] == {k: str(p) for k, p in paths.items()},
        "Assessment submission differs from original conversion/reviewer/request")
    require(sub["worker"] == record(assessment_gate.__file__)
        and sub["driver"] == record(assessment_gate.assessment.__file__)
        and sub["driver"]["sha256"] == assessment_gate.ASSESSMENT_SHA,
        "Bound assessment worker/driver changed")
    check(sub["batch"])
    return sub, conversion, reviewer, request, run, paths, conversion_paths


def completed_assessment(fields, sub):
    require(fields.get("JobIDRaw") == str(sub["job_id"]) and fields.get("State") == "COMPLETED"
        and fields.get("ExitCode") == "0:0" and fields.get("NodeList") == "bizon"
        and fields.get("AllocCPUS") == "8"
        and fields.get("ReqMem") in {"64G", "64Gn", "65536M", "65536Mn"},
        "Assessment is not successfully completed in its bound resource envelope")


def completed_gate(gate, sub, execution_ref):
    require(gate.get("schema") == "conversion_gated_native_qfo_assessment_v1"
        and gate.get("status") == "native_qfo_assessment_process_succeeded_pending_independent_admission"
        and gate.get("source") == sub["worker"] and gate.get("driver") == sub["driver"]
        and gate.get("conversion_submission") == sub["conversion_submission"]
        and gate.get("job_id") == str(sub["job_id"])
        and gate.get("native_job_id") == sub["native_job_id"]
        and gate.get("reviewer_job_id") == sub["reviewer_job_id"]
        and gate.get("conversion_job_id") == sub["conversion_job_id"]
        and gate.get("index") == sub["index"] and gate.get("cell") == sub["cell"]
        and gate.get("output_namespaces") == sub["output_namespaces"]
        and gate.get("execution") == execution_ref
        and gate.get("assessment_driver_invoked") is True and gate.get("endpoint_process_completed") is True
        and all(gate.get(k) is False for k in ("accuracy_admitted", "native_inference_reexecuted",
            "automatic_retry", "next_identity_authorized", "publication_ready")),
        "Require the successful bound assessment gate, not merely endpoint files")


def output_paths(run, native_job):
    return dict(gate=ROOT / "benchmarks/work" / f"native_factorial_qfo_admission_gate_{native_job}",
                admission=ROOT / "benchmarks/results/full_native_qfo_admission_v1" / run["cell"])


def execute(submission_ref, worker_sha256):
    source = record(__file__)
    require(source["sha256"] == worker_sha256, "Admission gate worker changed")
    sub, conversion, reviewer, request, run, assessment_paths, conversion_paths = submission_binding(submission_ref)
    paths = output_paths(run, sub["native_job_id"])
    for path in paths.values():
        require(path.is_absolute() and path.resolve() == path, "Require direct admission namespaces")
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    require(os.environ.get("SLURM_CPUS_PER_TASK") == "2" and os.environ.get("SLURM_JOB_ID", "").isdigit()
        and int(os.environ["SLURM_JOB_ID"]) not in
            {sub[k] for k in ("job_id", "native_job_id", "reviewer_job_id", "conversion_job_id")},
        "Require a separate scheduled two-CPU admission job")
    directory = paths["gate"]
    directory.mkdir(parents=True, exist_ok=False)
    report = dict(schema="assessment_gated_native_qfo_admission_v1", status="validating_after_assessment",
        source=source, assessment_submission=submission_ref, job_id=os.environ["SLURM_JOB_ID"],
        native_job_id=sub["native_job_id"], reviewer_job_id=sub["reviewer_job_id"],
        conversion_job_id=sub["conversion_job_id"], assessment_job_id=sub["job_id"],
        index=run["index"], cell=run["cell"], output_namespaces={k: str(p) for k, p in paths.items()},
        validator_invoked=False, accuracy_admitted=False, native_inference_reexecuted=False,
        automatic_retry=False, next_identity_authorized=False, publication_ready=False,
        started_monotonic_ns=time.monotonic_ns())
    try:
        text, scheduler = admission.accounting(sub["job_id"], include_memory=True)
        completed_assessment(scheduler, sub)
        report.update(assessment_accounting=text, assessment_scheduler=scheduler,
                      runtime=assessment_gate.conversion_gate.runtime_environment(reviewer))
        gate_ref = record(assessment_paths["gate"] / "results.json")
        execution_ref = record(assessment_paths["cwd"] / "results.json")
        gate = read(gate_ref)
        completed_gate(gate, sub, execution_ref)
        conversion_gate_ref, pairs_ref = [record(p / "results.json") for p in conversion_paths]
        assessment_gate.completed_gate(read(conversion_gate_ref), conversion, pairs_ref)
        require(gate["conversion_gate"] == conversion_gate_ref and gate["pairs"] == pairs_ref,
                "Assessment gate differs from actual bound conversion output")
        validator_ref = record(admission.__file__)
        require(validator_ref["sha256"] == ADMISSION_SHA, "Original independent validator changed")
        context = assessment_gate.conversion_gate.capacity()
        require(context["available_memory_bytes"] >= 32 * 2**30
            and context["available_disk_bytes"] >= 32 * 2**30, "Unsafe admission memory/disk capacity")
        refs = [source, submission_ref, gate_ref, execution_ref, conversion_gate_ref, pairs_ref, validator_ref,
            sub["worker"], sub["driver"], sub["batch"], sub["conversion_submission"],
            sub["reviewer_submission"], sub["request"], sub["plan"]]
        for ref in refs:
            check(ref)
        report.update(status="invoking_unchanged_independent_validator", validator_invoked=True,
            assessment_gate=gate_ref, execution=execution_ref, conversion_gate=conversion_gate_ref,
            pairs=pairs_ref, validator=validator_ref, resource_context=context)
        save(directory / "preflight.json", report)
        admitted = admission.admit(ROOT, pairs_ref, str(sub["conversion_job_id"]), str(sub["job_id"]), paths["admission"])
        admitted_ref = record(paths["admission"] / "results.json")
        report["original_admission"] = admitted_ref
        require(read(admitted_ref) == admitted
            and admitted.get("schema") == "full_native_factorial_qfo_admission_v1"
            and admitted.get("status") == "full_native_factorial_qfo_assessment_admitted"
            and admitted.get("source") == validator_ref and admitted.get("pairs_manifest") == pairs_ref
            and admitted.get("execution_report") == execution_ref and admitted.get("accuracy_admitted") is True
            and admitted.get("native_index") == run["index"] and admitted.get("cell") == run["cell"]
            and admitted.get("native_job_id") == sub["native_job_id"]
            and admitted.get("conversion") == read(pairs_ref)
            and admitted.get("participant") == "ohmm_qfo_full_native_" + run["cell"]
            and all(admitted.get(k) is False for k in ("publication_ready", "automatic_retry", "next_identity_authorized"))
            and all(admitted.get("scheduler", {}).get(k) == scheduler[k] for k in
                ("JobIDRaw", "State", "ExitCode", "NodeList", "AllocCPUS")),
            "Unexpected original independent admission identity or scope")
        for ref in refs:
            check(ref)
        report.update(status="assessment_gated_native_qfo_accuracy_admitted", accuracy_admitted=True)
    except BaseException as error:
        report.update(status="assessment_gated_native_qfo_admission_failed_retained",
                      error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_monotonic_ns"] = time.monotonic_ns()
        report["limitations"] = [
            "Original independent validator rechecks scientific/runtime/output/trace/FAS membership and arithmetic; no substitute validator.",
            "Before score export require this wrapper's successful terminal accounting and final accuracy-admitted gate, not a residual original report alone.",
            "Future output digests are observed only after producer completion and checked against bound submissions and gates.",
            "Point-score admission is not paired uncertainty, independent biological validation or publication readiness.",
            "Shared-host postprocessing outside inference timing; safe capacity does not establish isolation or whole-job budget adequacy.",
            "No inference, endpoint execution, automatic retry, timing repair or next-native-identity authorization."]
        save(directory / "results.json", report)
    return record(directory / "results.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assessment-submission", type=Path, required=True)
    parser.add_argument("--assessment-submission-sha256", required=True)
    parser.add_argument("--worker-sha256", required=True)
    args = parser.parse_args()
    ref = record(args.assessment_submission)
    require(ref["sha256"] == args.assessment_submission_sha256, "Assessment submission checksum differs")
    print(json.dumps(execute(ref, args.worker_sha256), sort_keys=True))

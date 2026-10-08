"""Allocation-aware successor to the retained preflight-failed runtime component."""

import argparse
import json
from pathlib import Path
import time

from benchmark_tools import review_native11_postterminal_runtime as kernel
from benchmark_tools.native_factorial_allocated_execution import amendment, validate_request, verify_terminal
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.review_allocated_native_factorial_attempt import bind_session
from benchmark_tools.review_native_factorial_attempt import scheduler_fields
from benchmark_tools.run_native_factorial_cost import ROOT, read
from benchmark_tools.validate_native_factorial_outputs import Evidence, require


SCHEMA = "native11_postterminal_runtime_component_v2"
KERNEL_SHA = "5c2e9c80c2c1be1f2972a914c08dcd4669909f755d73b445756f5d3144e6d8e0"
DESTINATION = ROOT / "benchmarks/work/native11_postterminal_runtime_review_20261008_v2"


def terminal_gate(terminal, request_ref):
    fields = scheduler_fields(terminal)
    require(fields.get("JobState", fields.get("State")) == "COMPLETED"
        and fields.get("ExitCode") == "0:0", "Require successful native outcome")
    if terminal["source"] == "live_controller":
        require(fields.get("Comment") == request_ref["sha256"], "Native comment differs")


def review(destination, expected_source_sha256):
    destination = Path(destination)
    require(destination == DESTINATION and destination.is_absolute()
        and destination.resolve() == destination and not destination.exists(),
        "Require fresh fixed v2 component destination")
    source, kernel_ref = record(__file__), record(kernel.__file__)
    require(source["sha256"] == expected_source_sha256
        and kernel_ref["sha256"] == KERNEL_SHA, "Prospective component sources changed")
    request_ref, classification_ref = record(kernel.REQUEST), record(kernel.CLASSIFICATION)
    require(request_ref["sha256"] == kernel.REQUEST_SHA
        and classification_ref["sha256"] == kernel.CLASSIFICATION_SHA, "Original bindings differ")
    destination.mkdir(parents=True, exist_ok=False)
    evidence = Evidence()
    try:
        classification = read(classification_ref)
        require(classification.get("schema") == "native11_postterminal_runtime_addition_classification_v1"
            and classification.get("classification_only") is True
            and classification.get("private_runtime_comparison_equal") is True
            and all(classification.get(k) is False for k in ("current_full_runtime_revalidated",
                "original_ordinary_full_review_success", "full_review_admitted", "next_identity_authorized")),
            "Require bound nonadmitting classification")
        for ref in [classification_ref, *classification["checked_records"],
                    classification["prior_inventory_source"]]:
            check(ref)
            evidence.bind(ref["path"])
        times = kernel.chronology(classification, read(classification["failed_review_observation"]))
        request = read(request_ref)
        execution, plan = amendment(request["amendment"])
        validate_request(request, request["amendment"], execution, 23985)
        run = plan["runs"][request["index"]]
        require(run["index"] == 11 and run["cell"] == "p1_c1_r0", "Wrong native identity")
        # The historical verifier expects a different job name, not this allocated route.
        terminal = verify_terminal(23985)
        terminal_gate(terminal, request_ref)
        save(destination / "scheduler.json", terminal)
        scheduler_ref = record(destination / "scheduler.json")
        for path, digest in kernel.SOURCES.items():
            ref = record(ROOT / path)
            require(ref["sha256"] == digest, "Frozen runtime kernel changed")
            evidence.bind(ref["path"])
        for ref in [source, kernel_ref, request_ref, request["plan"], request["amendment"],
                    plan["baseline"], *execution["new_sources"], *plan["helper_sources"], *plan["evidence"]]:
            check(ref)
            evidence.bind(ref["path"])
        session = Path(plan["panel_root"]) / "sessions" / "run_11"
        result = evidence.json(session / "result.json")
        verification = evidence.json(Path(run["output_root"]) / "verification.json")
        bind_session(result, verification, request_ref, request, request["amendment"], run)
        prior_ref = classification["original_runtime"]
        prior = read(prior_ref)
        evidence.bind(prior_ref["path"])
        report = kernel.replay_runtime(plan, run, session, verification, prior, evidence,
            classification["result"]["additional_entries"], times)
        report.update(schema=SCHEMA, source=source, runtime_kernel_source=kernel_ref,
            request=request_ref, classification=classification_ref, prior_runtime=prior_ref,
            plan=request["plan"], amendment=request["amendment"], scheduler=scheduler_ref,
            index=11, job_id=23985, cell=run["cell"], observed_unix_ns=time.time_ns(),
            terminal_verifier_source=record(ROOT / "benchmark_tools/native_factorial_allocated_execution.py"))
        report["evidence"] = evidence.finish()
        check(source)
        check(kernel_ref)
        check(scheduler_ref)
        save(destination / "runtime.json", report)
        return record(destination / "runtime.json")
    except Exception as error:
        save(destination / "failure.json", dict(schema=SCHEMA, status="runtime_component_failed",
            error_type=type(error).__name__, error=str(error), source=source,
            runtime_kernel_source=kernel_ref, request=request_ref, classification=classification_ref,
            full_review_admitted=False, next_identity_authorized=False,
            automatic_retry=False, publication_ready=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sha256", required=True)
    args = parser.parse_args()
    print(json.dumps(review(DESTINATION, args.source_sha256), sort_keys=True))

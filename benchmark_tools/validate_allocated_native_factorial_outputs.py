"""New production entry gate, unchanged native scientific semantic validation."""

import argparse
import json
from pathlib import Path

from benchmark_tools.native_factorial_allocated_execution import (
    ROOT, SCOPE, amendment, validate_request, verify_terminal, native_command)
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_native_factorial_cost import read
from benchmark_tools.run_allocated_native_factorial_cost import native_placement
from benchmark_tools.validate_native_factorial_outputs import require, validate_semantics


def validate(request_ref):
    request = read(request_ref)
    amendment_ref = request["amendment"]
    execution, plan = amendment(amendment_ref)
    plan_ref = execution["historical_plan"]
    validate_request(request, amendment_ref, execution, request["job_id"])
    terminal = verify_terminal(request["job_id"])
    fields = terminal["verified"].get("fields", terminal["verified"])
    require(fields.get("JobState", fields.get("State")) == "COMPLETED" and fields["ExitCode"] == "0:0",
            "Require actual successful allocated-native completion")
    if terminal["source"] == "live_controller":
        require(fields.get("Comment") == request_ref["sha256"], "Terminal allocated request comment differs")
    run = plan["runs"][request["index"]]
    baseline = read(plan["baseline"])
    root = Path(run["output_root"])
    execution_ref = record(root / "native_execution.json")
    native = read(execution_ref)
    ready_ref = record(root / "measurement/ready.json")
    ready = read(ready_ref)
    allowed = native_placement(ready, native["placement"], request["job_id"], native["parent_pid"])
    require(native.get("schema") == "allocated_native_factorial_execution_v1"
        and native.get("status") == "native_factorial_completed_pending_output_review"
        and native.get("plan") == plan_ref and native.get("amendment") == amendment_ref
        and native.get("index") == run["index"] and native.get("cell") == run["cell"]
        and native.get("allocated_ready") == ready_ref and native.get("native_cpu_ids") == allowed
        and native.get("native_order") == run["native_order"] and native.get("automatic_retry") is False
        and native.get("source") == record(ROOT / "benchmark_tools/run_allocated_native_factorial_cost.py"),
        "Native allocated execution provenance differs")
    evidence = [request_ref, amendment_ref, plan_ref, plan["baseline"], execution_ref, ready_ref,
        *execution["new_sources"], *plan["helper_sources"], *run["inputs"]]
    for ref in evidence:
        check(ref)
    context = dict(run, input_directory=str(root / "input"), cpu=32, threads_per_worker=4,
        command=native_command(amendment_ref, run, baseline, metrics=True), cwd=baseline["core_root"],
        aligner=baseline["tool_entrypoints"]["mafft"]["absolute_path"],
        tree_builder=baseline["tool_entrypoints"]["FastTree"]["absolute_path"])
    result = validate_semantics(context)
    preparation_ref = record(root / "preparation.json")
    preparation = read(preparation_ref)
    require(preparation["status"] == "fresh_factorial_inputs_prepared"
        and preparation["gene_ownership_sha256"] == result["gene_ownership_sha256"]
        and preparation["per_species_counts"] == result["per_species_counts"]
        and preparation["genes"] == run["genes"] and native["factors"] == result["factors"],
        "Prepared inputs or frozen native factors differ from scientific output")
    evidence.append(preparation_ref)
    for ref in evidence:
        check(ref)
    return dict(result, schema="allocated_native_factorial_output_review_v1",
        semantic_validator_source=result["source"], source=record(__file__),
        index=run["index"], job_id=request["job_id"], plan=plan_ref, amendment=amendment_ref,
        request=request_ref, native_cpu_ids=allowed, allocated_ready=ready_ref, scheduler=terminal,
        execution_scope=SCOPE, evidence=evidence, terminal_reviewed=False,
        terminal_scheduler_confirmed=True, uncontended_timing=False,
        contention_distortion="unknown_potentially_tool_dependent")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--request-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output review already exists")
    ref = record(args.request)
    require(ref["sha256"] == args.request_sha256, "Request checksum differs")
    report = validate(ref)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")

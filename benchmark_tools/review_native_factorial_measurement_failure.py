"""Review completed scientific outputs after a diagnosed cadence failure; never retry."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_native_factorial_cadence_failure import ERROR
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.review_native_factorial_attempt import bind_session, runtime_review, save, scheduler_fields
from benchmark_tools.run_native_factorial_cost import ROOT, SCOPE, read, validate_plan, validate_request, verify_terminal
from benchmark_tools.validate_native_factorial_outputs import Evidence, require, validate_semantics


def failure_scope(audit, request_ref, request, run, terminal, wrapper, done, environment):
    fields = scheduler_fields(terminal)
    require(fields.get("JobState", fields.get("State")) == "FAILED" and fields.get("ExitCode") == "1:0",
        "Keep the original failed allocation outcome")
    require(audit.get("schema") == "native_factorial_cadence_failure_audit_v1"
        and audit.get("status") == "retained_measurement_cadence_failure_reproduced"
        and audit.get("request") == request_ref and audit.get("plan") == request["plan"]
        and audit.get("job_id") == request["job_id"] and audit.get("index") == run["index"]
        and audit.get("cell") == run["cell"] and audit.get("primary_resources") is None
        and all(audit.get(k) is False for k in ("full_resource_replay", "scientific_timings_admitted",
            "native_outputs_validated", "next_identity_authorized", "original_receipts_rewritten",
            "inference_reexecuted", "automatic_retry")), "Require the original non-admitting cadence audit")
    census = audit["census"]
    failures = census.get("failures", [])
    require(failures and census.get("failure_counts") == {ERROR: len(failures)}
        and census.get("unchanged_cadence_bounds_s") == [.5, 1.5]
        and all(f.get("error") == ERROR and (f["wall_s"] < .5 or f["wall_s"] > 1.5)
                for f in failures), "Do not recover an unclassified or relaxed measurement failure")
    require(wrapper.get("status") == "verified_wrapper_failed" and wrapper.get("error") == ERROR
        and type(done.get("exit_code")) is int and done["exit_code"] == 0
        and done.get("timed_out") is False, "Require completed native command and failed measurement wrapper")
    pressure = environment.get("pressure_review", {})
    require(audit.get("retained_pressure_review") == pressure,
        "Cadence diagnosis and original pressure evidence differ")
    require(environment.get("job_id") == request["job_id"] and environment.get("index") == run["index"]
        and environment.get("execution_scope") == SCOPE and environment.get("uncontended_timing") is False
        and environment.get("sampled_process_policy_satisfied") is True
        and environment.get("sampled_environment_policy_satisfied") is False
        and environment.get("background_cpu_used_for_eligibility") is False
        and not environment.get("failures") and pressure.get("native_pressure_role") == "diagnostic_only"
        and pressure.get("pressure_thresholds_used_for_eligibility") is False
        and pressure.get("sampled_pressure_evidence_satisfied") is False
        and pressure.get("failures") == {"pressure_sample_period_exceeded": 1},
        "Require the retained classified cadence-only environment failure, not unsafe or unclassified evidence")


def reuse_runtime(plan, run, session, wrapper, runtime_ref, evidence):
    original = read(runtime_ref)
    require(original.get("status") == "runtime_brackets_and_lookup_replayed"
        and original.get("continuous_runtime_integrity_established") is False,
        "Require original successful runtime review")
    original_copy = dict(original)
    original_copy.pop("fresh_terminal_check_wall_s")
    replayed = runtime_review(plan, run, session, wrapper, evidence,
        tree_checker=lambda specifications: original["fresh_terminal_inventory"])
    replayed.pop("fresh_terminal_check_wall_s")
    require(replayed == original_copy, "Original runtime brackets, lookup or manifest bindings differ")
    evidence.bind(runtime_ref["path"])
    return dict(status="retained_runtime_brackets_bound_and_replayed", original_runtime_review=runtime_ref,
        current_manifest_and_lookup_bindings_checked=True, fresh_runtime_tree_recheck=False,
        original_terminal_inventory_reused=True, continuous_runtime_integrity_established=False,
        limitations=["Replay reuses the original terminal tree verdict; no new whole-runtime tree audit is claimed.",
            "Manifest/lookup/input bindings are checked; future launches retain their fresh runtime gates."])


def review(request_ref, audit_ref, failed_ref, destination):
    destination = Path(destination)
    request, audit, failed = read(request_ref), read(audit_ref), read(failed_ref)
    plan_ref = request["plan"]
    plan = read(plan_ref)
    run = validate_plan(plan)[request["index"]]
    validate_request(request, plan_ref, request["job_id"])
    root = Path(run["output_root"])
    session = Path(plan["panel_root"]) / "sessions" / f"run_{run['index']:02d}"
    require(destination.is_absolute() and destination.resolve() == destination
        and destination.is_relative_to(ROOT) and not destination.is_relative_to(root)
        and not destination.is_relative_to(session) and not destination.exists()
        and not destination.is_symlink(), "Require separate fresh failure-disposition directory")
    require(failed.get("status") == "terminal_factorial_review_failed"
        and failed.get("request") == request_ref and failed.get("plan") == plan_ref
        and failed.get("job_id") == request["job_id"] and failed.get("index") == run["index"]
        and failed.get("terminal_reviewed") is False and failed.get("next_identity_authorized") is False
        and "lineage_report.json" in failed.get("error", "")
        and failed.get("source") == record(ROOT / "benchmark_tools/review_native_factorial_attempt.py"),
        "Require the unchanged original failed reviewer")
    require(audit.get("evidence", [])[3] ==
        record(ROOT / "benchmark_tools/audit_native_factorial_cadence_failure.py"),
        "Require the unchanged actual cadence-audit source")
    evidence = Evidence()
    for ref in [request_ref, audit_ref, failed_ref, plan_ref, plan["baseline"],
                *audit["evidence"], *plan["helper_sources"], *plan["evidence"]]:
        require(record(evidence.bind(ref["path"])) == ref, "Bound failure evidence differs")
    source = record(__file__)
    terminal = verify_terminal(request["job_id"])
    result = evidence.json(session / "result.json")
    wrapper = evidence.json(root / "verification.json")
    bind_session(result, wrapper, request_ref, request, plan_ref, run)
    done = evidence.json(root / "measurement/done.json")
    environment_ref = result["environment_review"]
    require(environment_ref["path"] == str(root / "measurement/process_stream_review.json"),
        "Wrong retained environment report")
    environment = read(environment_ref)
    evidence.bind(environment_ref["path"])
    failure_scope(audit, request_ref, request, run, terminal, wrapper, done, environment)
    runtime_ref = record(Path(failed_ref["path"]).parent / "runtime.json")
    runtime = reuse_runtime(plan, run, session, wrapper, runtime_ref, evidence)
    baseline = read(plan["baseline"])
    command = [baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"],
        str(ROOT / "benchmark_tools/run_native_factorial_cost.py"), "--native", "--plan", plan_ref["path"],
        "--plan-sha256", plan_ref["sha256"], "--index", str(run["index"])]
    context = dict(run, input_directory=str(root / "input"), cpu=32, threads_per_worker=4,
        command=command, cwd=baseline["core_root"],
        aligner=baseline["tool_entrypoints"]["mafft"]["absolute_path"],
        tree_builder=baseline["tool_entrypoints"]["FastTree"]["absolute_path"])
    outputs = validate_semantics(context)
    native = evidence.json(root / "native_execution.json")
    metrics = evidence.json(root / "metrics.json")
    preparation = evidence.json(root / "preparation.json")
    require(native.get("status") == "native_factorial_completed_pending_output_review"
        and native.get("plan") == plan_ref and native.get("index") == run["index"]
        and native.get("cell") == run["cell"] and native.get("factors") == outputs["factors"]
        and native.get("counts") == metrics["counts"] and native.get("stages") == sorted(metrics["stages"])
        and native.get("native_order") == run["native_order"] and native.get("automatic_retry") is False,
        "Scientific terminal receipt or frozen factors differ")
    require(preparation.get("status") == "fresh_factorial_inputs_prepared"
        and preparation.get("genes") == run["genes"]
        and preparation.get("gene_ownership_sha256") == outputs["gene_ownership_sha256"]
        and preparation.get("per_species_counts") == outputs["per_species_counts"],
        "Recovered output universe differs from prepared input provenance")
    for ref in outputs["checked_files"]:
        require(record(evidence.bind(ref["path"])) == ref, "Recovered semantic evidence differs")
    evidence.bind(source["path"])
    checked = evidence.finish()
    destination.mkdir(parents=True, exist_ok=False)
    runtime_readback = save(destination / "runtime_readback.json", runtime)
    output_ref = save(destination / "outputs.json", outputs)
    fields = scheduler_fields(terminal)
    return save(destination / "review.json", dict(schema="native_factorial_measurement_failure_review_v1",
        status="native_scientific_outputs_recovered_measurement_failure_retained", job_id=request["job_id"],
        index=run["index"], cell=run["cell"], dataset=run["dataset"], repeat=run["repeat"],
        plan=plan_ref, request=request_ref, source=source, scheduler=terminal,
        scheduler_state=fields.get("JobState", fields.get("State")), scheduler_exit_code=fields["ExitCode"],
        original_failed_review=failed_ref, cadence_diagnosis=audit_ref, environment_report=environment_ref,
        runtime_readback=runtime_readback, outputs=output_ref, evidence=checked,
        terminal_reviewed=True, next_identity_authorized=True,
        next_identity_scope="Only the next different frozen identity, after its unchanged fresh launch gates; never this attempt again.",
        native_outputs_validated=True, accuracy_evaluated=False, primary_resources_replayed=False,
        shared_host_resources_reviewed=False, resources=None, scientific_timings_admitted=False,
        eligible_for_timing_comparison=False, scheduler_success=False, native_command_success=True,
        original_receipts_rewritten=False, inference_reexecuted=False, automatic_retry=False,
        execution_scope=SCOPE, uncontended_timing=False,
        contention_distortion="unknown_potentially_tool_dependent",
        limitations=["Distinct failed-measurement disposition, never a successful terminal resource review.",
            "Scientific outputs are semantically checked, not accuracy-scored or independently biologically validated.",
            "Resource evidence remains ineligible; missing reports are not fabricated and cadence criteria remain unchanged.",
            "Retained runtime brackets are reused with checked bindings, not a new whole-runtime or continuous audit.",
            "Next identity may proceed only with its original fresh capacity, accounting, runtime and input gates.",
            "Shared-host effects remain unknown; no corrected time, isolated ranking, retry or fastest-repeat selection."]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "audit", "failed-review", "output-directory"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("request-sha256", "audit-sha256", "failed-review-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    refs = [record(path) for path in (args.request, args.audit, args.failed_review)]
    require([r["sha256"] for r in refs] ==
        [args.request_sha256, args.audit_sha256, args.failed_review_sha256], "Explicit recovery binding differs")
    print(json.dumps(review(*refs, args.output_directory.resolve()), sort_keys=True))

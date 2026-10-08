"""Explicit conversion binding for the new composed review, not a schema shim."""

from pathlib import Path

from benchmark_tools import review_native11_composed_attempt as composed
from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting
from benchmark_tools.derive_threadripper_resources import SCOPES
from benchmark_tools.native_factorial_allocated_execution import amendment, validate_request, verify_terminal, SCOPE
from benchmark_tools.prepare_native_factorial_qfo_pairs import conversion_kind
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.review_native_factorial_attempt import scheduler_fields
from benchmark_tools.run_native_factorial_cost import ROOT, read
from benchmark_tools.validate_native_factorial_outputs import require


COMPOSED_SOURCE_SHA = "1fa148071e8d080ea081140f517e955c16db46f989574a44183fc58e1b366410"
HELD = ROOT / "benchmark_tools/results/native11_composed_review_held_24034_20261008_v1.json"
HELD_SHA = "d6785db02c31034873c26b6f8350dac2672beee2a4387b87b5035c7c4b7b1427"
RELEASE = ROOT / "benchmark_tools/results/native11_composed_review_release_24034_20261008_v1.json"
RELEASE_SHA = "64258ad0fc6fd3714e3d044cc0afedc9e52594f01a76a47efd6ec84b3ba692d8"


def admit_composed_review(review, request_ref, request, run, producer):
    require(run.get("index") == 11 and run.get("cell") == "p1_c1_r0"
        and run.get("dataset") == "qfo_corrected"
        and review.get("schema") == composed.SCHEMA and review.get("status") == "native_success"
        and review.get("request") == request_ref and review.get("plan") == request["plan"]
        and review.get("amendment") == request["amendment"] and request.get("job_id") == 23985
        and review.get("job_id") == 23985 and review.get("review_job_id") == 24034
        and review.get("index") == 11 and review.get("cell") == run["cell"]
        and review.get("dataset") == run["dataset"] and type(review.get("repeat")) is int
        and review["repeat"] == run["repeat"]
        and review.get("scheduler_state") == "COMPLETED" and review.get("scheduler_exit_code") == "0:0"
        and review.get("execution_scope") == SCOPE and review.get("resource_scopes") == SCOPES
        and all(review.get(k) is True for k in ("terminal_reviewed", "composed_full_review_complete",
            "native_outputs_validated", "primary_resources_replayed", "shared_host_resources_reviewed"))
        and all(review.get(k) is False for k in ("current_original_inventory_equality",
            "original_ordinary_full_review_success", "next_identity_authorized", "downstream_adoption_complete",
            "accuracy_evaluated", "scientific_timings_admitted", "uncontended_timing", "automatic_retry", "publication_ready"))
        and review.get("historical_failures_retained") == [23986, 24033]
        and set(review.get("reviews", {})) == {"runtime", "resources", "environment", "outputs_or_failure"},
        "Require explicit successful composed native11 review, not ordinary-schema substitution")
    require((producer.get("JobIDRaw"), producer.get("State"), producer.get("ExitCode"),
             producer.get("NodeList"), producer.get("AllocCPUS")) ==
            ("24034", "COMPLETED", "0:0", "bizon", "2")
        and producer.get("ReqMem") in {"128G", "128Gn", "131072M", "131072Mn"},
        "Composed producer is not successfully terminal in its declared envelope")
    return conversion_kind(run)


def native_binding(request_ref, review_ref):
    request, review = read(request_ref), read(review_ref)
    require(request_ref["sha256"] == composed.runtime_kernel.REQUEST_SHA
        and review_ref["path"] == str(composed.DESTINATION / "review.json"),
        "Require original request and separate composed result destination")
    execution, plan = amendment(request["amendment"])
    validate_request(request, request["amendment"], execution, 23985)
    run = plan["runs"][request["index"]]
    producer_raw, producer = accounting(24034, include_memory=True)
    kind = admit_composed_review(review, request_ref, request, run, producer)
    source_ref = record(composed.__file__)
    require(source_ref["sha256"] == COMPOSED_SOURCE_SHA and review.get("source") == source_ref,
        "Composed reviewer source differs")
    held_ref, release_ref = record(HELD), record(RELEASE)
    require(held_ref["sha256"] == HELD_SHA and release_ref["sha256"] == RELEASE_SHA,
        "Composed producer submission/release changed")
    held, release = read(held_ref), read(release_ref)
    require(held.get("job_id") == release.get("job_id") == 24034
        and held.get("held_comparison_passed") is True and held.get("submission_count") == 1
        and release.get("returncode") == 0 and release.get("release_count") == 1
        and release.get("held") == held_ref and held["references"]["worker"] == source_ref,
        "Composed producer release/source binding differs")
    allocation = read(review["allocation"])
    require(allocation["fields"]["JobId"] == "24034", "Composed allocation belongs to another job")
    composed.allocation_gate(allocation["stdout"], 24034, COMPOSED_SOURCE_SHA)
    terminal = verify_terminal(23985)
    composed.runtime_v2.terminal_gate(terminal, request_ref)
    fields = scheduler_fields(terminal)
    require(fields.get("JobState", fields.get("State")) == "COMPLETED", "Native no longer successful")
    runtime = read(review["reviews"]["runtime"])
    require(runtime.get("schema") == "native11_composed_runtime_gate_v1"
        and runtime.get("source") == source_ref and runtime.get("runtime_component") == record(composed.RUNTIME)
        and runtime.get("current_original_inventory_equality") is False
        and runtime.get("current_original_entries_unchanged") is True
        and runtime.get("current_private_inventory_equality") is True,
        "Composed runtime gate provenance differs")
    environment = read(review["reviews"]["environment"])
    resources = read(review["reviews"]["resources"])
    require(environment.get("sampled_environment_evidence_valid") is True
        and environment.get("amendment") == request["amendment"]
        and resources.get("primary") == review["resources"] and resources.get("primary_scopes") == SCOPES,
        "Composed resource/environment reports differ")
    output_ref = review["reviews"]["outputs_or_failure"]
    output = read(output_ref)
    require(output_ref == record(composed.OUTPUT) and output_ref["sha256"] == composed.OUTPUT_SHA
        and output.get("native_outputs_validated") is True and output.get("request") == request_ref
        and output.get("plan") == request["plan"] and output.get("amendment") == request["amendment"]
        and output.get("native_cpu_ids") == review["native_cpu_ids"]
        and output.get("allocated_ready") == record(Path(run["output_root"]) / "measurement/ready.json"),
        "Composed semantic output/placement differs")
    records = [request_ref, review_ref, source_ref, held_ref, release_ref, request["plan"], request["amendment"],
        review["allocation"], review["scheduler"], review["semantic_producer"], review["resource_replay"],
        *review["reviews"].values(), *review["evidence"], *output["checked_files"], *output["evidence"],
        output["source"], output["semantic_validator_source"], output["allocated_ready"],
        *held["references"].values(), plan["baseline"], *execution["new_sources"],
        *plan["helper_sources"], *plan["evidence"]]
    for ref in records:
        check(ref)
    return request, execution, plan, run, review, output, kind, terminal, records, dict(
        review_producer_accounting=producer_raw, review_producer_scheduler=producer,
        composed_schema_preserved=True, original_review_translated=False, next_identity_authorized=False)

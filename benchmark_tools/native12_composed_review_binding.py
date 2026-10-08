"""Bind a successful final native review without translating its schema."""

import ast
import os
from pathlib import Path
import re

from benchmark_tools import review_native12_composed_attempt as reviewer
from benchmark_tools import run_native12_composed_cost as executor
from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting
from benchmark_tools.derive_threadripper_resources import SCOPES
from benchmark_tools.native_factorial_allocated_execution import SCOPE
from benchmark_tools.prepare_native_factorial_qfo_pairs import conversion_kind
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.review_native_factorial_attempt import scheduler_fields
from benchmark_tools.run_native_factorial_cost import ROOT, read
from benchmark_tools.validate_native_factorial_outputs import require


REQUEST_SHA = "70bde16142da526de6f198aaf255a224b0f4e10c2bd79f98f0b1fd7d70269f89"
REVIEWER_SHA = "08b170ba1ee5aabd5665dd12fc359fc7f5ffe3d0834b0afb2e73289fa61ed649"
REVIEW_BATCH = ROOT / "benchmark_tools/results/native12_composed_terminal_review_20261008_v1.sh"
REVIEW_BATCH_SHA = "0439a2c230909e3d060bb2c3b72e7bab6862e1109521203241f8165b948d6e93"


def producer_envelope(raw, job):
    pairs = re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", raw)
    fields = dict(pairs)
    expected = dict(JobId=str(job), JobName="ohmm_native12_review", JobState="PENDING",
        Reason="JobHeldUser", Partition="gpu", ReqNodeList="bizon", NumCPUs="2",
        NumTasks="1", MinMemoryNode="128G", TimeLimit="06:00:00", Requeue="0", Restarts="0",
        Command=str(REVIEW_BATCH), WorkDir=str(ROOT), Comment=REVIEWER_SHA)
    expected["CPUs/Task"] = "2"
    require(type(job) is int and job > 24036 and len(fields) == len(pairs)
        and len([line for line in raw.splitlines() if line.strip()]) == 1
        and all(fields.get(key) == value for key, value in expected.items())
        and fields.get("NumNodes") in {"1", "1-1"}
        and fields.get("UserId", "").endswith("(" + str(os.getuid()) + ")")
        and fields.get("Dependency") in {None, "(null)"}
        and not any(key in fields for key in ("ArrayJobId", "ArrayTaskId", "HetJobId", "HetJobOffset")),
        "Final review producer held envelope differs")
    return fields


def admit_review(review, request_ref, request, run, producer, producer_job):
    require(run.get("index") == 12 and type(run["index"]) is int
        and run.get("cell") == "p1_c1_r1" and run.get("dataset") == "qfo_corrected"
        and request.get("schema") == executor.REQUEST_SCHEMA and request.get("job_id") == 24036
        and review.get("schema") == reviewer.SCHEMA and review.get("status") == "native_success"
        and review.get("request") == request_ref and review.get("plan") == request["plan"]
        and review.get("amendment") == request["amendment"]
        and review.get("job_id") == request["job_id"] and type(review.get("index")) is int
        and review["index"] == 12 and review.get("cell") == run["cell"]
        and review.get("dataset") == run["dataset"] and type(review.get("repeat")) is int
        and review["repeat"] == run["repeat"]
        and review.get("scheduler_state") == "COMPLETED" and review.get("scheduler_exit_code") == "0:0"
        and review.get("execution_scope") == SCOPE and review.get("resource_scopes") == SCOPES
        and all(review.get(key) is True for key in ("terminal_reviewed", "native_outputs_validated",
            "primary_resources_replayed", "shared_host_resources_reviewed", "prospective_current_inventory_equality"))
        and all(review.get(key) is False for key in ("original_review_translated",
            "current_original_os_inventory_equality", "continuous_runtime_integrity_established",
            "next_identity_authorized", "accuracy_evaluated", "scientific_timings_admitted",
            "uncontended_timing", "automatic_retry", "publication_ready"))
        and review.get("historical_review_failures_retained") == [23986, 24033]
        and set(review.get("reviews", {})) == {"runtime", "resources", "environment", "outputs_or_failure"},
        "Require successful explicit final-identity review, not another schema or cell")
    require(type(producer_job) is int and producer_job > 24036
        and (producer.get("JobIDRaw"), producer.get("State"), producer.get("ExitCode"),
             producer.get("NodeList"), producer.get("AllocCPUS")) ==
            (str(producer_job), "COMPLETED", "0:0", "bizon", "2")
        and producer.get("ReqMem") in {"128G", "128Gn", "131072M", "131072Mn"},
        "Final review producer is not successfully terminal in its declared envelope")
    return conversion_kind(run)


def native_binding(request_ref, review_ref, producer_job, held_ref, release_ref):
    require(request_ref["path"] == str(executor.REQUEST) and request_ref["sha256"] == REQUEST_SHA
        and review_ref["path"] == str(reviewer.DESTINATION / "review.json"),
        "Require actual final request and separate final review destination")
    request = read(request_ref)
    request, context, history = executor.execution_binding(request_ref, request["job_id"])
    execution, plan = context[1:3]
    executor.held_gate(request["held_scheduler"]["stdout"], request["job_id"])
    run, review = plan["runs"][12], read(review_ref)
    raw, producer = accounting(producer_job, include_memory=True)
    kind = admit_review(review, request_ref, request, run, producer, producer_job)
    source_ref, batch_ref = record(reviewer.__file__), record(REVIEW_BATCH)
    require(source_ref["sha256"] == REVIEWER_SHA and review.get("source") == source_ref
        and batch_ref["sha256"] == REVIEW_BATCH_SHA
        and source_ref in request["new_sources"] and batch_ref in request["new_sources"],
        "Final reviewer source or batch differs from inference-bound sources")
    held, release = read(held_ref), read(release_ref)
    require(held.get("schema") == "native12_composed_review_submission_v1"
        and release.get("schema") == "native12_composed_review_release_v1"
        and held.get("job_id") == release.get("job_id") == producer_job
        and held.get("submission_count") == 1 and held.get("held_comparison_passed") is True
        and release.get("returncode") == 0 and release.get("release_count") == 1
        and release.get("held") == held_ref
        and held.get("destination") == str(reviewer.DESTINATION)
        and held.get("native_inference_reexecuted") is False and held.get("automatic_retry") is False
        and held.get("references", {}).get("request") == request_ref
        and held["references"].get("worker") == source_ref and held["references"].get("batch") == batch_ref,
        "Final review submission/release/source provenance differs")
    producer_fields = producer_envelope(held["controller"]["stdout"], producer_job)
    producer_logs = [record(Path(producer_fields[key])) for key in ("StdOut", "StdErr")]
    require(ast.literal_eval(Path(producer_logs[0]["path"]).read_text().strip()) == review_ref,
        "Final review producer stdout does not bind this completed review")
    require(review["scheduler"]["path"] == str(reviewer.DESTINATION / "scheduler.json")
        and review["resource_replay"]["path"] == str(reviewer.DESTINATION / "resource_replay_summary.json")
        and all(review["reviews"][key]["path"] == str(reviewer.DESTINATION / filename)
            for key, filename in (("runtime", "runtime.json"), ("resources", "resources.json"),
                ("environment", "environment.json"), ("outputs_or_failure", "outputs_or_failure.json"))),
        "Final review component destination differs")
    terminal = reviewer.verify_terminal(request["job_id"])
    fields = scheduler_fields(terminal)
    require(fields.get("JobState", fields.get("State")) == "COMPLETED"
        and fields.get("ExitCode") == "0:0", "Final inference no longer successfully terminal")
    if terminal["source"] == "live_controller":
        require(fields.get("Comment") == request_ref["sha256"], "Final inference request comment differs")
    runtime = read(review["reviews"]["runtime"])
    require(runtime.get("schema") == "native12_composed_runtime_review_v1"
        and runtime.get("status") == "fresh_runtime_brackets_and_lookup_replayed"
        and runtime.get("runtime_basis") == request["runtime_basis"]
        and runtime.get("current_original_os_inventory_equality") is False
        and runtime.get("prospective_current_inventory_equality") is True
        and runtime.get("continuous_runtime_integrity_established") is False
        and set(runtime.get("phases", {})) == {"before", "after"},
        "Final runtime review provenance differs")
    resources = read(review["reviews"]["resources"])
    environment = read(review["reviews"]["environment"])
    replayed = read(review["resource_replay"])
    require(resources.get("primary") == review["resources"] and resources.get("primary_scopes") == SCOPES
        and environment.get("sampled_environment_evidence_valid") is True
        and environment.get("amendment") == request["amendment"]
        and all(environment.get(key) is False for key in ("uncontended_timing",
            "background_cpu_used_for_eligibility", "pressure_thresholds_used_for_eligibility"))
        and replayed.get("schema") == "native12_composed_resource_replay_summary_v1"
        and replayed.get("full_replay_executed") is True
        and replayed.get("measured_matches_retained_wrapper") is True
        and replayed.get("source") == record(ROOT / "benchmark_tools/replay_allocated_threadripper_scaling.py"),
        "Final environment or full resource replay binding differs")
    output = read(review["reviews"]["outputs_or_failure"])
    require(output.get("schema") == "native12_composed_output_review_v1"
        and output.get("source") == source_ref and output.get("native_outputs_validated") is True
        and output.get("accuracy_evaluated") is False
        and output.get("index") == 12 and output.get("job_id") == request["job_id"]
        and output.get("cell") == run["cell"] and output.get("request") == request_ref
        and output.get("plan") == request["plan"] and output.get("amendment") == request["amendment"]
        and output.get("native_cpu_ids") == review["native_cpu_ids"]
        and output.get("allocated_ready") == record(Path(run["output_root"]) / "measurement/ready.json")
        and output.get("semantic_validator_source") == record(ROOT / "benchmark_tools/validate_native_factorial_outputs.py")
        and type(output.get("phylogeny", {}).get("native_pair_rows")) is int
        and output["phylogeny"]["native_pair_rows"] >= 0,
        "Final semantic output/placement/native pair binding differs")
    records = [request_ref, review_ref, source_ref, batch_ref, held_ref, release_ref, request["plan"],
        request["amendment"], review["scheduler"], review["resource_replay"], *review["reviews"].values(),
        *review["evidence"], *output["checked_files"], *output["evidence"], output["source"],
        output["semantic_validator_source"], output["allocated_ready"], *held["references"].values(),
        *replayed["evidence"], plan["baseline"], *execution["new_sources"], *request["new_sources"],
        *plan["helper_sources"], *plan["evidence"], *history["evidence"], *producer_logs]
    for ref in records:
        check(ref)
    return request, execution, plan, run, review, output, kind, terminal, records, dict(
        review_producer_job_id=producer_job, review_producer_accounting=raw,
        review_producer_scheduler=producer, review_held=held_ref, review_release=release_ref,
        composed_schema_preserved=True, original_review_translated=False, next_identity_authorized=False)

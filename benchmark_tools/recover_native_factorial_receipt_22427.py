"""Recover scientific outputs of 22427, never rewrite its failed native outcome."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.collect_factorial_resources import FIXED_INPUTS
from benchmark_tools.link_factorial_scaling_resources import partition, compare
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.review_native_factorial_attempt import save, scheduler_fields
from benchmark_tools.run_native_factorial_cost import ROOT, read, validate_plan, verify_terminal
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.validate_native_factorial_outputs import Evidence, require, validate_semantics


BROKEN_EXECUTOR_SHA = "07a926e38ba49144b56bf9a7d23478d3dd050fc3a14b98e1ce8889d1f98a8d20"
CREATE_ONCE_WRITER_SHA = "c00eae1c44c86f0ee02c72e4634cf4e5cc70e2656a8ebfc0f9070769f6f7f97b"


def receipt_failure(initial, pending, metrics, run, plan_ref, log):
    require(initial.get("status") == "native_factorial_running"
        and pending.get("status") == "native_factorial_completed_pending_output_review"
        and metrics.get("status") == "complete" and pending.get("counts") == metrics["counts"]
        and pending.get("stages") == sorted(metrics["stages"]), "Incomplete scientific completion evidence")
    require(all(pending.get(k) == v for k, v in initial.items() if k != "status")
        and initial.get("plan") == plan_ref and initial.get("index") == run["index"]
        and initial.get("cell") == run["cell"] and initial.get("native_order") == run["native_order"]
        and initial.get("automatic_retry") is False, "Pending receipt changes execution identity")
    path = Path(run["output_root"]) / "native_execution.json"
    error = f"FileExistsError: [Errno 17] File exists: '{path}.pending' -> '{path}'"
    require(log.rstrip().splitlines()[-1] == error and log.count("FileExistsError:") == 1
        and "save(report, receipt)" in log, "Failure is not the identified completion-receipt collision")


def score_orthobench(run, outputs, evidence):
    relative, digest = FIXED_INPUTS["ob_results"]
    source = evidence.bind(ROOT / relative)
    require(evidence.files[str(source)]["sha256"] == digest, "Frozen OrthoBench reference/score manifest differs")
    prior = json.loads(source.read_text())
    references, uncertain = {}, {}
    for ref in prior["references"]:
        path = evidence.bind(ref["path"])
        require(evidence.files[str(path)] == ref, "Frozen OrthoBench reference checksum differs")
        target = uncertain if path.parent.name == "low_certainty_assignments" else references
        require(path.name not in target, "Duplicate RefOG evidence")
        target[path.name] = set(path.read_text().split())
    require(sorted(references) == [f"RefOG{i:03d}.txt" for i in range(1, 71)], "Require the frozen seventy RefOGs")
    prediction = Path(run["output_root"]) / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
    predicted = partition(evidence.bind(prediction), "space_separated_groups")
    require(len(predicted[0]) == outputs["input_genes"], "Scored output universe differs")
    old_ref = prior["predictions"][run["cell"]]
    old_path = evidence.bind(old_ref["path"])
    require(evidence.files[str(old_path)] == old_ref, "Retained original prediction checksum differs")
    old = partition(old_path, "space_separated_groups")
    old_score = score_partition(old[1], references, uncertain)
    for key in ("f_score", "precision", "recall"):
        require(abs(old_score[key] - prior["point_estimates_percent"][run["cell"]][key]) < 1e-10,
                "Audited reference statistic does not reproduce")
    score = score_partition(predicted[1], references, uncertain)
    reference_genes = set().union(*references.values())
    delta = compare(old, predicted, reference_genes)
    scorer = evidence.bind(ROOT / "benchmark_tools/score_orthobench_partition.py")
    return dict(schema="recovered_native_orthobench_score_v1", dataset="OrthoBench", cell=run["cell"],
        score=score, score_fraction={k: score[k] / 100 for k in ("f_score", "precision", "recall")},
        original_score_percent=old_score, original_partition_comparison=delta,
        metric="official equal-RefOG weighted pair statistic with low-certainty exclusions",
        development_exposed=True, independent_validation=False, timing_success_established=False,
        reference_manifest=evidence.files[str(source)], scorer=evidence.files[str(scorer)])


def recover(request_ref, review_ref, destination):
    request, terminal_review = read(request_ref), read(review_ref)
    require(request["job_id"] == terminal_review["job_id"] == 22427 and request["index"] == 0
        and terminal_review["request"] == request_ref and terminal_review["plan"] == request["plan"]
        and terminal_review["terminal_reviewed"] is True
        and terminal_review["primary_resources_replayed"] is True
        and terminal_review["shared_host_resources_reviewed"] is True
        and terminal_review["status"] == "native_failure_retained", "Require the fully reviewed 22427 failure")
    terminal = verify_terminal(22427)
    fields = scheduler_fields(terminal)
    require(fields.get("JobState", fields.get("State")) == "FAILED" and fields["ExitCode"] == "1:0",
            "Do not relabel the terminal scheduler failure")
    plan_ref = request["plan"]
    plan = read(plan_ref)
    run = validate_plan(plan)[0]
    baseline = read(plan["baseline"])
    require(run["cell"] == "p0_c0_r0" and run["dataset"] == "orthobench", "Wrong recovery identity")
    evidence = Evidence()
    for name, digest in (("run_native_factorial_cost.py", BROKEN_EXECUTOR_SHA),
                         ("probe_dgx_step_separation.py", CREATE_ONCE_WRITER_SHA)):
        path = evidence.bind(ROOT / "benchmark_tools" / name)
        require(evidence.files[str(path)]["sha256"] == digest, "Failure source differs from the diagnosed revision")
    for ref in [request_ref, review_ref, terminal_review["source"], terminal_review["resource_replay"],
                *terminal_review["evidence"], *terminal_review["reviews"].values()]:
        require(evidence.files.get(ref["path"], ref) == ref, "Review evidence binding differs")
        path = evidence.bind(ref["path"])
        require(evidence.files[str(path)] == ref, "Retained review evidence changed")
    root = Path(run["output_root"])
    initial = evidence.json(root / "native_execution.json")
    pending = evidence.json(root / "native_execution.json.pending")
    metrics = evidence.json(root / "metrics.json")
    log = evidence.bind(root / "measurement/native.log").read_text()
    receipt_failure(initial, pending, metrics, run, plan_ref, log)
    command = [baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"],
        str(ROOT / "benchmark_tools/run_native_factorial_cost.py"), "--native", "--plan", plan_ref["path"],
        "--plan-sha256", plan_ref["sha256"], "--index", "0"]
    context = dict(run, input_directory=str(root / "input"), cpu=32, threads_per_worker=4,
        command=command, cwd=baseline["core_root"],
        aligner=baseline["tool_entrypoints"]["mafft"]["absolute_path"],
        tree_builder=baseline["tool_entrypoints"]["FastTree"]["absolute_path"])
    outputs = validate_semantics(context)
    require(pending["factors"] == outputs["factors"], "Scientific factor binding differs")
    preparation = evidence.json(root / "preparation.json")
    require(preparation["gene_ownership_sha256"] == outputs["gene_ownership_sha256"]
        and preparation["per_species_counts"] == outputs["per_species_counts"], "Preparation provenance differs")
    for ref in outputs["checked_files"]:
        path = evidence.bind(ref["path"])
        require(evidence.files[str(path)] == ref, "Recovered semantic evidence changed")
    score = score_orthobench(run, outputs, evidence)
    destination = Path(destination)
    require(destination.is_absolute() and destination.resolve() == destination
        and destination.is_relative_to(ROOT) and not destination.is_relative_to(root), "Require separate recovery destination")
    destination.mkdir(parents=True, exist_ok=False)
    outputs_ref = save(destination / "outputs.json", outputs)
    score_ref = save(destination / "score.json", score)
    source = evidence.bind(__file__)
    return save(destination / "recovery.json", dict(schema="native_factorial_receipt_failure_recovery_v1",
        status="scientific_outputs_recovered_from_failed_wrapper", job_id=22427, index=0,
        cell=run["cell"], plan=plan_ref, request=request_ref, terminal_review=review_ref, scheduler=terminal,
        outputs=outputs_ref, score=score_ref, evidence=evidence.finish(),
        source=evidence.files[str(source)], native_outputs_validated=True, accuracy_evaluated=True,
        scheduler_success=False, native_command_success=False, timing_success_established=False,
        original_receipts_rewritten=False, inference_reexecuted=False, automatic_retry=False,
        resources=terminal_review["resources"], resource_scopes=terminal_review["resource_scopes"],
        timing_scope="failed native-wrapper attempt after complete scientific pipeline; not a clean-success timing",
        limitations=["The scheduler and native command remain failed; recovery does not erase their outcomes.",
            "Scientific outputs are independently checked and rescored without repeating inference.",
            "Development-exposed OrthoBench evidence, not independent generalization or superiority.",
            "Shared-host resource observations retain their unknown, potentially tool-dependent contention effects."]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "terminal-review", "output-directory"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("request-sha256", "terminal-review-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    request_ref, review_ref = record(args.request), record(args.terminal_review)
    require(request_ref["sha256"] == args.request_sha256 and review_ref["sha256"] == args.terminal_review_sha256,
            "Recovery evidence checksum differs")
    print(json.dumps(recover(request_ref, review_ref, args.output_directory.resolve()), sort_keys=True))

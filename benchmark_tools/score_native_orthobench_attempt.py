"""Score independently terminal-reviewed native OrthoBench factorial outputs.

No inference, timing correction, method selection or next-job release occurs.
Corrected QfO requires its separate pair conversion/native endpoint workflow.
"""

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.collect_factorial_resources import FIXED_INPUTS
from benchmark_tools.derive_threadripper_resources import SCOPES
from benchmark_tools.link_factorial_scaling_resources import compare, partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import ROOT, SCOPE, read, validate_plan, validate_request, verify_terminal
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.validate_native_factorial_outputs import Evidence, require


def admit_score(review, request_ref, request, run):
    require(run["dataset"] == "orthobench", "QfO needs the separate native endpoint assessment, not RefOG scoring")
    require(review.get("schema") == "native_factorial_terminal_review_v1"
        and review.get("request") == request_ref and review.get("plan") == request["plan"]
        and review.get("job_id") == request["job_id"] and review.get("index") == run["index"]
        and review.get("cell") == run["cell"] and review.get("dataset") == run["dataset"]
        and review.get("repeat") == run["repeat"] and review.get("status") == "native_success"
        and review.get("scheduler_state") == "COMPLETED" and review.get("scheduler_exit_code") == "0:0"
        and all(review.get(k) is True for k in ("terminal_reviewed", "primary_resources_replayed",
            "shared_host_resources_reviewed", "native_outputs_validated"))
        and review.get("execution_scope") == SCOPE and review.get("resource_scopes") == SCOPES
        and review.get("automatic_retry") is False and review.get("uncontended_timing") is False
        and set(review.get("reviews", {})) == {"runtime", "resources", "environment", "outputs_or_failure"},
        "Require the bound successful terminal/runtime/resource/environment/output review")


def scored_partitions(current, original, references, uncertain, old_point):
    require(current[0] == original[0], "Prediction universe differs from frozen original")
    reference_genes = set().union(*references.values())
    require(reference_genes <= current[0], "RefOG genes missing from native inference universe")
    old_score = score_partition(original[1], references, uncertain)
    for key in ("f_score", "precision", "recall"):
        value = old_point.get(key)
        require(type(value) in (int, float) and math.isfinite(value)
            and math.isclose(old_score[key], value, rel_tol=0, abs_tol=1e-10),
            "Frozen original weighted statistic does not reproduce: " + key)
    score = score_partition(current[1], references, uncertain)
    delta = compare(original, current, reference_genes)
    delta["native_only_groups"] = delta.pop("scaling_only_groups")
    return dict(score_percent=score,
        score_fraction={key: score[key] / 100 for key in ("f_score", "precision", "recall")},
        original_score_percent=old_score,
        difference_from_original_percentage_points={key: score[key] - old_score[key]
            for key in ("f_score", "precision", "recall")},
        canonical_partition_comparison=delta, reference_genes=len(reference_genes),
        covered_reference_genes=len(reference_genes & current[0]), inference_genes=len(current[0]))


def score_frozen(run, output_review, evidence):
    name, digest = FIXED_INPUTS["ob_results"]
    manifest_ref = record(ROOT / name)
    require(manifest_ref["sha256"] == digest, "Frozen OrthoBench reference/results manifest changed")
    manifest = read(manifest_ref)
    evidence.bind(manifest_ref["path"])
    references, uncertain = {}, {}
    for ref in manifest["references"]:
        check(ref)
        path = evidence.bind(ref["path"])
        target = uncertain if path.parent.name == "low_certainty_assignments" else references
        require(path.name not in target, "Duplicate frozen RefOG reference")
        target[path.name] = set(path.read_text().split())
    require(sorted(references) == [f"RefOG{i:03d}.txt" for i in range(1, 71)], "Require all seventy frozen RefOGs")
    reconciliation = run["cell"].endswith("r1")
    format = "root_hogs" if reconciliation else "space_separated_groups"
    suffix = "native/orthohmm_phylogeny/orthohmm_root_hogs.tsv" if reconciliation else "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
    predicted_path = Path(run["output_root"]) / suffix
    checked = {r["path"]: r for r in output_review["checked_files"]}
    require(str(predicted_path) in checked, "Chosen prediction was not independently output-validated")
    predicted_ref = checked[str(predicted_path)]
    check(predicted_ref)
    current = partition(evidence.bind(predicted_ref["path"]), format)
    original_ref = manifest["predictions"][run["cell"]]
    check(original_ref)
    original = partition(evidence.bind(original_ref["path"]), format)
    require(len(current[0]) == output_review["input_genes"] == run["genes"]
        and len(current[1]) == output_review["orthogroups"], "Scored universe/group count differs from admitted output")
    result = scored_partitions(current, original, references, uncertain,
        manifest["point_estimates_percent"][run["cell"]])
    return dict(result, reference_manifest=manifest_ref, prediction=predicted_ref,
        original_prediction=original_ref, prediction_format=format,
        prediction_semantics="root-HOG co-membership" if reconciliation else "final cluster co-membership",
        endpoint="official equal-RefOG weighted pair statistic with low-certainty exclusions",
        reference_families=70, development_exposed=True, independent_validation=False)


def score_attempt(request_ref, review_ref, destination):
    source_ref = record(__file__)
    request, review = read(request_ref), read(review_ref)
    plan = read(request["plan"])
    validate_request(request, request["plan"], request["job_id"])
    run = validate_plan(plan)[request["index"]]
    admit_score(review, request_ref, request, run)
    destination = Path(destination)
    require(destination.is_absolute() and destination.resolve() == destination
        and destination.is_relative_to(ROOT) and not destination.is_relative_to(Path(plan["panel_root"])),
        "Require separate direct score destination")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    terminal = verify_terminal(request["job_id"])
    fields = terminal["verified"].get("fields", terminal["verified"])
    require(fields.get("JobState", fields.get("State")) == "COMPLETED" and fields.get("ExitCode") == "0:0",
            "Fresh scheduler observation is not successful terminal")
    if terminal["source"] == "live_controller":
        require(fields.get("Comment") == request_ref["sha256"], "Current scheduler request comment differs")
    evidence = Evidence()
    pins = [request_ref, request["plan"], review_ref, review["source"], review["scheduler"],
        *review["reviews"].values(), plan["baseline"], *plan["helper_sources"], *run["inputs"]]
    require(review["source"] == record(ROOT / "benchmark_tools/review_native_factorial_attempt.py"),
            "Terminal reviewer source differs")
    for ref in pins:
        check(ref)
        evidence.bind(ref["path"])
    outputs = read(review["reviews"]["outputs_or_failure"])
    require(outputs.get("native_outputs_validated") is True and outputs.get("request") == request_ref
        and outputs.get("plan") == request["plan"] and outputs.get("job_id") == request["job_id"]
        and outputs.get("index") == run["index"] and outputs.get("cell") == run["cell"],
        "Independent output review identity differs")
    for ref in [*outputs["checked_files"], *outputs["evidence"]]:
        check(ref)
        evidence.bind(ref["path"])
    result = score_frozen(run, outputs, evidence)
    require(record(evidence.bind(__file__)) == source_ref, "Scorer source changed during evaluation")
    for name in ("score_native_orthobench_attempt.py", "score_orthobench_partition.py", "link_factorial_scaling_resources.py"):
        evidence.bind(ROOT / "benchmark_tools" / name)
    result.update(schema="native_factorial_orthobench_score_v1", status="terminal_native_orthobench_scored",
        job_id=request["job_id"], index=run["index"], cell=run["cell"], dataset="OrthoBench", repeat=run["repeat"],
        request=request_ref, plan=request["plan"], terminal_review=review_ref, scheduler=terminal,
        source=source_ref, evidence=evidence.finish(), native_outputs_validated=True, accuracy_evaluated=True,
        native_inference_reexecuted=False, automatic_retry=False, next_identity_authorized=False,
        publication_ready=False, resource_observation=review["resources"], resource_scopes=SCOPES,
        timing_scope=SCOPE, uncontended_timing=False, timing_disclosure=review["timing_disclosure"],
        limitations=["Development-exposed accuracy, not independent confirmation or selection-adjusted superiority.",
            "Root HOGs/cluster membership are scored here; native resolved ortholog pairs require a different endpoint.",
            "Comparison with cached outputs tests reproducibility, not a new method's biological contribution.",
            "Only bound terminal-reviewed successes use this route; separately recovered failures retain their failed outcomes.",
            "Shared-host time/resource observations may be distorted by unknown tool-dependent contention.",
            "No inference rerun, bootstrap claim, timing correction or automatic next-job release."])
    destination.mkdir(parents=True, exist_ok=False)
    path = destination / "score.json"
    save(path, result)
    return record(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "terminal-review", "output-directory"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("request-sha256", "terminal-review-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    request_ref, review_ref = record(args.request), record(args.terminal_review)
    require(request_ref["sha256"] == args.request_sha256 and review_ref["sha256"] == args.terminal_review_sha256,
            "Scoring request/review checksum differs")
    print(json.dumps(score_attempt(request_ref, review_ref, args.output_directory.absolute()), sort_keys=True))

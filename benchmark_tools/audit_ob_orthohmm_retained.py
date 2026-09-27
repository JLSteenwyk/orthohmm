"""Audit historical OrthoHMM OrthoBench records without substituting new runs."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

from benchmark_tools.run_qfo_order_replay import record, save
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.prepare_qfo_canonical_pairs import check_records


def resolve_records(base, items):
    result = []
    seen = set()
    for item in items:
        path = Path(item["path"])
        if path.is_absolute() or ".." in path.parts or str(path) in seen:
            raise ValueError("Invalid or duplicate relative manifest path")
        seen.add(str(path))
        result.append({**item, "path": str((base / path).resolve())})
    return result


def timing(metrics):
    for key in ("wall_s", "user_cpu_s", "system_cpu_s", "peak_process_tree_rss_bytes"):
        value = metrics[key]
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            raise ValueError("Invalid historical resource value")
    elapsed = metrics["finished_at_epoch_s"] - metrics["started_at_epoch_s"]
    if not math.isclose(elapsed, metrics["wall_s"], rel_tol=0, abs_tol=1e-5):
        raise ValueError("Historical elapsed time disagrees")
    if metrics["rss_measurement"] != "sampled_sum_of_linux_proc_tree_rss":
        raise ValueError("Unknown memory semantics")
    return {k: metrics[k] for k in ("wall_s", "user_cpu_s", "system_cpu_s",
                                   "peak_process_tree_rss_bytes", "rss_measurement")}


def source_evidence(repo, commit, items):
    rows = []
    for item in items:
        payload = subprocess.run(["git", "show", f"{commit}:{item['path']}"], cwd=repo,
                                 capture_output=True, check=False)
        digest = hashlib.sha256(payload.stdout).hexdigest() if payload.returncode == 0 else None
        rows.append(dict(recorded=item, commit=commit, blob_sha256=digest,
                         blob_matches=digest == item["sha256"] and len(payload.stdout) == item["bytes"]))
    return rows


def audit(repo, output):
    if output.exists():
        raise FileExistsError(output)
    results = repo / "benchmarks/results"
    high_path = results / "production_high_sensitivity_replay_profiles_seed4_committed_20260828/result.json"
    phylo_path = results / "production_phylogeny_satellite_v2_orthobench_20260902/metrics.json"
    high = read_frozen(high_path, "642cb044d01771070f06a6557ca5a0c0bde82f9ea68eb5797346d5fec1e9306f")
    phylo = read_frozen(phylo_path, "dbde87690a592a6bfa1adfbbb560106d593d5a1b5ae59febea73f34d37603f49")
    scores_path = repo / "benchmark_tools/results/current_benchmark_scores_20260926_v2/manifest.json"
    scores = read_frozen(scores_path, "8d56dae74721c4530ab2b86c5008b2d0ebac07878e574fd85b1ad84e4b908f07")
    score_rows = {r["key"]: r for r in scores["rows"]}
    harness = phylo["harness"]
    if phylo["status"] != "complete" or harness["exit_code"] != 0:
        raise ValueError("Incomplete historical inference")
    inputs = resolve_records(Path(phylo["metadata"]["fasta_directory"]), harness["input_manifest"])
    outputs = resolve_records(phylo_path.parent, harness["output_manifest"])
    if len(inputs) != 12:
        raise ValueError("Expected twelve input proteomes")
    high_outputs = [s["output"] for s in high["stages"]]
    high_prediction = score_rows["orthohmm_high_sensitivity"]["orthobench_retained_evidence"]["prediction_provenance"]
    phylo_prediction = score_rows["orthohmm_phylogeny_satellite_v2"]["orthobench_retained_evidence"]["prediction_provenance"]
    selected = [s for s in high["stages"] if s["label"] == "strict_profiles_refined"]
    if len(selected) != 1 or selected[0]["output"] != high_prediction or phylo_prediction not in outputs:
        raise ValueError("Metrics do not bind scored predictions")
    checked = [record(p) for p in (high_path, phylo_path, scores_path, Path(__file__))]
    checked += [high["input"], *high_outputs, *inputs, *outputs]
    check_records(checked)
    high_source = {**high["source"], "path": str(Path(high["source"]["path"]).relative_to(repo))}
    sources = source_evidence(repo, high["git"]["commit"], [high_source])
    sources += source_evidence(repo, harness["git_commit"], harness["source_manifest"])
    report = dict(status="retained_orthohmm_ob_provenance_audited", checked_records=checked,
        high_sensitivity=dict(command=high["command"], prediction=high_prediction, hit_input=high["input"],
            parameters=high["parameters"], counts=high["counts"], wall_s=high["wall_s"],
            peak_process_rss_gib=high["peak_process_rss_gib"], timings=high["timings"],
            scope="Cached-hit downstream replay, excluding initial search", input_fasta_history_proven=False),
        phylogeny=dict(command=phylo["command"], prediction=phylo_prediction,
            inputs=inputs, outputs=outputs, resources=timing(phylo), stages=phylo["stages"],
            metadata=phylo["metadata"], counts=phylo["counts"],
            scope="Historical full inference, shared host; scoring excluded"),
        historical_source_blobs=sources, all_recorded_sources_match_commits=all(r["blob_matches"] for r in sources),
        dirty_historical_worktrees=dict(high=high["git"]["dirty"], phylogeny=harness["git_dirty"]),
        publication_ready=False, controlled_comparative_resources=False, complete_transitive_provenance=False,
        limitations=["Retained metrics are not immutable proof of historical execution or absence of contention",
            "Cached replay and full inference timings are not comparable scopes",
            "Sampled process-tree RSS double-counts shared pages and is not physical unique memory",
            "Historical dirty checkouts and external runtime identities require separate audit",
            "No historical score or resource value is filled from fresh recovery"])
    check_records(checked)
    save(output, report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    audit(args.repo.resolve(), args.output.resolve())

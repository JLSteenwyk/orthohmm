"""Run prespecified complete-proteome OrthoBench profile-evidence pilots.

The historical weighting interface remains the default. Alignment and
selective HMMER select their separately prespecified experiments.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from Bio import SeqIO
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_orthobench_factorial import prepare_partition
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
from benchmark_tools.run_canonical_downstream import commands
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import (
    PRODUCTS, compare_golden, score_orthobench, scoring_records, validate_partition, write_json,
)
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import (
    candidate_control_matches, sources, timed_stage, validate_ownership,
)
from benchmark_tools.validate_profile_runtime import require_profile_runtime
from benchmark_tools.native_pair_turnover import native_pair_turnover
from orthohmm.accuracy import load_accuracy_checkpoint
from orthohmm.orthohmm import _expand_phylogeny_candidates


ARMS = ("none", "henikoff", "henikoff_cluster", "species_henikoff_cluster")
PROTOCOL = "hmm_evidence_campaign_20261018.json"
RETAINED = "benchmarks/results/canonical_downstream_r1_v1_20261014/results.json"
UNIVERSE = (251378, 12)


def experiment_settings(experiment):
    if experiment == "weighting":
        return ARMS, PROTOCOL
    if experiment == "alignment":
        return ("off", "strict_span_guard"), "hmm_aligned_profile_guard_20261018.json"
    if experiment == "selective_hmmer":
        return ("off", "builtin_reselect", "hmmer_score_only", "hmmer_reselect"), "hmm_selective_hmmer_20261018.json"
    if experiment == "projection_recovery":
        return ("legacy_off", "projection_off", "projection_builtin_reselect",
                "projection_hmmer_score_only", "projection_hmmer_reselect"), "hmm_projection_integrity_recovery_20261018.json"
    raise ValueError("Unknown profile-evidence experiment")


def recovery_configuration(mode):
    arms, _ = experiment_settings("projection_recovery")
    if mode not in arms:
        raise ValueError("Unknown projection recovery arm")
    if mode == "legacy_off":
        return "legacy", "off"
    return "residue_preserving", mode.removeprefix("projection_")


def advancement(control, candidate):
    """Apply the frozen point gate, without claiming statistical equivalence."""
    fields = ("f1", "precision", "recall")
    for row in (control, candidate):
        if set(row) != set(fields) or any(
            isinstance(row[key], bool) or not isinstance(row[key], (int, float))
            or not math.isfinite(row[key]) or not 0 <= row[key] <= 1 for key in fields
        ):
            raise ValueError("Require finite normalized F1/precision/recall")
    delta = {key: candidate[key] - control[key] for key in fields}
    return {"eligible": delta["f1"] > 0
            and candidate["precision"] >= control["precision"] - .005
            and candidate["recall"] >= control["recall"] - .005,
            "deltas": delta, "scope": "Development point gate, not equivalence or independent confirmation"}


def selected_baseline(repo):
    retained_path = repo / RETAINED
    retained = json.loads(retained_path.read_text())
    candidates = [row for row in retained["arms"]
                  if row["dataset"] == "OrthoBench" and row["label"] == "canonical"]
    if len(candidates) != 1 or candidates[0]["status"] != "independent_content_readback_verified":
        raise ValueError("Require exactly one admitted canonical OrthoBench control")
    row = candidates[0]
    check(retained["execution_plan"])
    original = json.loads(Path(retained["execution_plan"]["path"]).read_text())
    case = next(case for case in original["cases"] if case["dataset"] == "OrthoBench")
    native = Path(original["output"]) / "OrthoBench/canonical/native"
    checkpoint_pin = row["ordering"]["manifest"]
    check(checkpoint_pin)
    baseline = {
        "dataset": "OrthoBench", "fasta_directory": case["fasta_directory"],
        "baseline_output": str(native), "checkpoint": str(Path(checkpoint_pin["path"]).parent),
        "checkpoint_manifest": checkpoint_pin, "fastas": case["fastas"],
        "candidate": row["candidate"]["candidate_partition"],
        "constraints": row["candidate"]["membership_constraints"],
        "aligner": case["aligner"], "tree_builder": case["tree_builder"],
        "golden": {Path(pin["path"]).name: pin for pin in row["outputs"]},
        "expected": {"f1": row["scores"]["OrthoBench F1"],
                     "precision": row["details"]["precision"], "recall": row["details"]["recall"]},
        "evidence": [record(retained_path), retained["execution_plan"]],
    }
    baseline["golden"]["species_tree.rooted.nwk"] = record(native / "orthohmm_phylogeny/species_tree.rooted.nwk")
    return baseline


def prepare(repo, source, output, experiment="weighting"):
    if output.exists() or source == repo:
        raise ValueError("Require fresh output and a frozen scientific checkout")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
    arms, protocol_name = experiment_settings(experiment)
    protocol = source / "benchmark_tools" / protocol_name
    spec = json.loads(protocol.read_text())
    declared = (tuple(arm["label"] for arm in spec["first_experiment"]["arms"])
                if experiment == "weighting" else tuple(spec["recovery_panel"]["ordered_arms"])
                if experiment == "projection_recovery" else tuple(spec["arms"]))
    if declared != arms:
        raise ValueError("Unexpected prespecified profile-evidence arms")
    case = selected_baseline(repo)
    if len(case["fastas"]) != 12:
        raise ValueError("Require all12 OrthoBench proteomes")
    inputs = [record(protocol), *case["evidence"], *case["fastas"], case["candidate"],
              case["constraints"], *case["golden"].values(), *scoring_records(repo, "OrthoBench"),
              *[record(path) for path in sorted(Path(case["checkpoint"]).iterdir()) if path.is_file()]]
    annotation_protocol = source / "benchmark_tools/hmm_stockholm_annotation_recovery_20261019.json"
    if experiment == "projection_recovery" and annotation_protocol.is_file():
        inputs.append(record(annotation_protocol))
    phylo = Path(case["baseline_output"]) / "orthohmm_phylogeny"
    for directory in ("alignments", "gene_trees", "checkpoints", "species_tree_inference"):
        inputs.extend(record(path) for path in sorted((phylo / directory).rglob("*")) if path.is_file())
    inputs.extend(record(shutil.which(case[key]) or case[key]) for key in ("aligner", "tree_builder"))
    hmmer_tools = {}
    if experiment in ("selective_hmmer", "projection_recovery"):
        for tool in ("hmmbuild", "hmmsearch"):
            binary = Path(shutil.which(tool) or tool).resolve(strict=True)
            help_text = subprocess.check_output([str(binary), "-h"], text=True, timeout=10)
            if not any(line.startswith("# HMMER 3.4 ") for line in help_text.splitlines()):
                raise ValueError("Require explicit HMMER3.4 pilot tools")
            hmmer_tools[tool] = record(binary)
        inputs.extend(hmmer_tools.values())
    inputs = list({pin["path"]: pin for pin in inputs}.values())
    source_pins = sources(source) + [record(source / "benchmark_tools/run_profile_weighting_pilot.slurm")]
    for pin in source_pins:
        relative = Path(pin["path"]).relative_to(source)
        if relative.suffix != ".so":
            blob = subprocess.check_output(["git", "show", revision + ":" + str(relative)], cwd=source)
            if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
                raise ValueError("Frozen source differs from committed Git bytes: " + str(relative))
    original_source = repo / "benchmarks/work/full_search_hillclimb_r1_v3_20261012"
    for name in ("pair_align.so", "hmm_viterbi.so", "kmer_prefilter.so"):
        relative = Path("orthohmm/search/csrc") / name
        if record(source / relative)["sha256"] != record(original_source / relative)["sha256"]:
            raise ValueError("Require unchanged original CPU scoring libraries")
    scientific_runtime = runtime()
    for pin in (*inputs, *source_pins, *scientific_runtime["files"]):
        check(pin)
    plan = {"status": "prepared_unrun", "source_root": str(source), "source_revision": revision,
            "repo": str(repo), "output": str(output), "cpu": 16, "arms": list(arms),
            "experiment": experiment,
            "case": case, "inputs": inputs, "sources": source_pins, "runtime": scientific_runtime,
            "profile_runtime": require_profile_runtime(source), "tool_path": os.environ["PATH"],
            "scope": "Complete-proteome downstream HMM/profile inference from immutable canonical initial-search hits, with content-compatible tree-cache reuse. Shared-host incremental timing, not a fresh-search or isolated full-pipeline speed comparison.",
            "defaults_promoted": False}
    if experiment in ("selective_hmmer", "projection_recovery"):
        plan["hmmer_tools"] = hmmer_tools
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", plan)
    return plan


def arm_commands(source, plan, folder, mode):
    experiment = plan.get("experiment", "weighting")
    arms, _ = experiment_settings(experiment)
    if mode not in arms:
        raise ValueError("Unexpected weighting arm")
    case = plan["case"]
    graph, native = commands(source, plan, case, folder, Path(case["checkpoint"]), case["checkpoint_manifest"]["sha256"])
    if experiment == "projection_recovery":
        policy, selective = recovery_configuration(mode)
        if policy != "legacy":
            graph.extend(["--profile-projection-policy", policy])
        if selective != "off":
            graph.extend(["--selective-profile-mode", selective])
    elif mode != arms[0]:
        flag = {"weighting": "--profile-sequence-weighting", "alignment": "--profile-alignment-policy",
                "selective_hmmer": "--selective-profile-mode"}[experiment]
        graph.extend([flag, mode])
    return graph, native


def validate_alignment_evidence(iteration, strict_edges):
    evidence, sources = iteration.get("alignment_evidence", []), iteration.get("alignment_profiles", [])
    if not evidence or not sources:
        raise ValueError("Complete-proteome guard requires nonempty actual alignment evidence")
    by_profile = {row["cluster_id"]: row["member_ids"] for row in sources}
    if len(by_profile) != len(sources) or any(not 3 <= len(members) <= 200 for members in by_profile.values()):
        raise ValueError("Missing/duplicate source-profile member evidence")
    seen, edges, rejected = set(), {}, 0
    for row in evidence:
        gene, anchor, profile = row["gene_id"], row["anchor_gene_id"], row["profile_cluster_id"]
        if gene in seen or gene == anchor or profile not in by_profile or row["anchor_name"] not in by_profile[profile]:
            raise ValueError("Alignment evidence does not describe unique selected source anchors")
        seen.add(gene)
        m, t, matched = row["model_length"], row["target_length"], row["matched_known_residues"]
        if (not all(type(row[field]) is int for field in (
                "gene_id", "anchor_gene_id", "profile_cluster_id", "model_length", "target_length",
                "score", "model_start", "model_end", "target_start", "target_end", "matched_known_residues"))
                or gene < 0 or anchor < 0 or m <= 0 or t <= 0 or row["score"] <= 0
                or not 0 <= row["model_start"] < row["model_end"] <= m
                or not 0 <= row["target_start"] < row["target_end"] <= t
                or not 0 <= matched <= min(row["model_end"] - row["model_start"], row["target_end"] - row["target_start"])
                or not math.isfinite(row["anchor_score"]) or row["anchor_score"] <= 0):
            raise ValueError("Invalid alignment evidence coordinates, counts or scores")
        model_covered, target_covered = 5 * matched >= 4 * m, 5 * matched >= 4 * t
        keep = model_covered or target_covered
        category = ("both_sides_covered" if model_covered and target_covered else
                    "asymmetric_coverage" if keep else "partial_shared_region")
        expected = {"aligned_model_coverage": matched / m, "aligned_target_coverage": matched / t,
                    "bounding_model_coverage": (row["model_end"] - row["model_start"]) / m,
                    "bounding_target_coverage": (row["target_end"] - row["target_start"]) / t}
        if (type(row["passes_guard"]) is not bool or row["passes_guard"] != keep
                or row["category"] != category or row["coordinate_system"] != "zero_based_half_open"
                or any(row[key] != value for key, value in expected.items())):
            raise ValueError("Alignment classification or coverage does not reproduce from counts")
        if keep:
            key = min(gene, anchor), max(gene, anchor)
            edges[key] = max(edges.get(key, -math.inf), row["anchor_score"])
        else:
            rejected += 1
    if len(edges) != strict_edges or set(by_profile) != {row["profile_cluster_id"] for row in evidence}:
        raise ValueError("Guarded unique-edge/source-profile counts disagree with evidence")
    return {"alignment_traced": len(evidence), "alignment_rejected": rejected,
            "alignment_retained": len(evidence) - rejected}


def validate_activation(graph, mode, experiment="weighting", *, folder=None, gene_names=None, hmmer_tools=None,
                        projection_policy="legacy"):
    if experiment == "projection_recovery":
        policy, selective = recovery_configuration(mode)
        actual = validate_activation(graph, selective, "selective_hmmer", folder=folder,
                                     gene_names=gene_names, hmmer_tools=hmmer_tools, projection_policy=policy)
        return {**actual, "arm": mode, "projection_policy": policy}
    arms, _ = experiment_settings(experiment)
    if mode not in arms:
        raise ValueError("Unknown profile-evidence arm")
    parameters, counts = graph["parameters"], graph["counts"]
    expected_weighting = mode if experiment == "weighting" else "none"
    expected_alignment = mode if experiment == "alignment" else "off"
    if (parameters.get("profile_sequence_weighting") != expected_weighting
            or parameters.get("profile_alignment_policy", "off") != expected_alignment
            or parameters.get("profile_projection_policy", "legacy") != projection_policy
            or parameters.get("profile_iterations") != 1
            or parameters.get("profile_expansion") is not True
            or parameters.get("matrix") != "BLOSUM62"
            or parameters.get("cpm_resolution") != .1
            or parameters.get("leiden_seed") != 4
            or parameters.get("jackknife_profile_thresholds") is not False
            or parameters.get("jackknife_single_copy_profiles") is not False
            or parameters.get("profile_min_species") != 1
            or parameters.get("singleton_tie_policy") != "first_hit"
            or counts.get("profiles_built", 0) <= 0):
        raise ValueError("Graph/profile intervention differs from frozen settings or is inactive")
    iterations = graph.get("profile_iterations", [])
    if len(iterations) != 1 or iterations[0]["profiles_built"] != counts["profiles_built"]:
        raise ValueError("Missing profile iteration evidence")
    effective = iterations[0]["profile_effective_counts"]
    if expected_weighting == "none":
        if effective:
            raise ValueError("Unweighted control unexpectedly records weighted profiles")
    elif (len(effective) != counts["profiles_built"]
          or len({item["cluster_id"] for item in effective}) != len(effective)
          or any(not math.isfinite(item["effective_sequence_count"])
                 or not 0 < item["effective_sequence_count"] <= 200 + 1e-8 for item in effective)):
        raise ValueError("Invalid per-profile weighting activation/effective counts")
    result = {"mode": mode, "profiles_built": counts["profiles_built"],
              "weighted_profiles": len(effective), "strict_profile_edges": counts["strict_profile_edges"]}
    if expected_alignment == "strict_span_guard":
        result.update(validate_alignment_evidence(iterations[0], counts["strict_profile_edges"]))
    elif iterations[0].get("alignment_evidence") or iterations[0].get("alignment_profiles"):
        raise ValueError("Off control unexpectedly has alignment intervention evidence")
    if experiment == "selective_hmmer":
        from benchmark_tools.validate_selective_profile_evidence import validate
        result.update(validate(graph, mode, folder, gene_names, hmmer_tools))
    elif graph.get("selective_profile_evidence") or parameters.get("selective_profile_mode", "off") != "off":
        raise ValueError("Unexpected selective profile intervention in another experiment")
    return result


def normalized(row):
    return {"f1": row["scores"]["OrthoBench F1"],
            "precision": row["details"]["precision"], "recall": row["details"]["recall"]}


def run(plan_path):
    plan = json.loads(plan_path.read_text())
    source, repo, output = (Path(plan[key]) for key in ("source_root", "repo", "output"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_profile_weighting_pilot.py":
        raise ValueError("Execute frozen pilot source")
    experiment = plan.get("experiment", "weighting")
    arms, _ = experiment_settings(experiment)
    if tuple(plan["arms"]) != arms or plan["cpu"] != 16 or (output / "results.json").exists():
        raise ValueError("Unexpected arm/budget or existing attempt; preserve previous runs")
    checked = [record(plan_path), *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in checked:
        check(pin)
    if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Changed runtime or CUDA visibility")
    case = plan["case"]
    owners = {}
    for pin in case["fastas"]:
        for sequence in SeqIO.parse(pin["path"], "fasta"):
            if sequence.id in owners or not sequence.seq:
                raise ValueError("Duplicate or empty input gene")
            owners[sequence.id] = Path(pin["path"]).name
    if (len(owners), len(set(owners.values()))) != UNIVERSE:
        raise ValueError("Require complete OrthoBench gene/species universe")
    checkpoint = Path(case["checkpoint"])
    checkpoint_audit = audit(checkpoint, case["checkpoint_manifest"]["sha256"])
    names, species, q, t, scores = load_accuracy_checkpoint(checkpoint, verify=True)
    ownership = validate_ownership(names, species, owners)
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "PYTHONHASHSEED": "0",
           "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
           "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "16"}
    report = {"status": "running", "execution_plan": record(plan_path), "job_id": os.environ.get("SLURM_JOB_ID"),
              "source_revision": plan["source_revision"], "arms": [], "scope": plan["scope"],
              "experiment": experiment,
              "checkpoint_audit": checkpoint_audit, "ownership": ownership,
              "host_before": host_activity(), "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", report)
    selective_experiment = experiment in ("selective_hmmer", "projection_recovery")
    recovery = experiment == "projection_recovery"
    control, projection_control, measurement_control, candidate_digest, control_stages = None, None, False, None, None
    for mode in arms:
        selective_mode = recovery_configuration(mode)[1] if recovery else mode
        row = {"dataset": "OrthoBench", "label": mode, "status": "running", "scores": {}}
        report["arms"].append(row)
        if mode != arms[0] and control is None:
            row["status"] = "not_run_failed_control_gate"
            write_json(output / "results.json", report)
            continue
        if recovery and mode not in ("legacy_off", "projection_off") and projection_control is None:
            row["status"] = "not_run_failed_projection_control_gate"
            write_json(output / "results.json", report)
            continue
        if selective_experiment and selective_mode == "hmmer_reselect" and not measurement_control:
            row["status"] = "not_run_failed_measurement_control_gate"
            write_json(output / "results.json", report)
            continue
        folder = output / mode
        folder.mkdir(parents=True, exist_ok=False)
        started = time.monotonic()
        write_json(output / "results.json", report)
        try:
            for pin in checked:
                check(pin)
            graph, native_command = arm_commands(source, plan, folder, mode)
            write_json(folder / "preflight.json", {"graph": graph, "native": native_command})
            row["graph_execution"] = timed_stage(graph, folder, "graph", source, env)
            if row["graph_execution"]["returncode"]:
                raise RuntimeError("Graph/profile replay failed")
            graph_report = json.loads((folder / "graph_metrics.json").read_text())
            activation_options = ({"folder": folder, "gene_names": names, "hmmer_tools": plan["hmmer_tools"]}
                                  if selective_experiment else {})
            row["activation"] = validate_activation(graph_report, mode, experiment, **activation_options)
            stage_signature = [(stage["label"], stage["output"]["bytes"], stage["output"]["sha256"])
                               for stage in graph_report["stages"]] if selective_experiment else None
            if selective_experiment:
                if selective_mode == "off":
                    control_stages = stage_signature
                else:
                    digest = row["activation"]["candidate_digest"]
                    if candidate_digest is not None and candidate_digest != digest:
                        raise ValueError("Candidate generation differs across active arms")
                    candidate_digest = digest
                    if selective_mode == "hmmer_score_only" and stage_signature != control_stages:
                        raise ValueError("HMMER measurement-only graph products differ from off")
            row["graph_report"] = record(folder / "graph_metrics.json")
            for stage in graph_report["stages"]:
                check(stage["output"])
                read_partition(Path(stage["output"]["path"]), set(names))
            before, candidate_started = resource.getrusage(resource.RUSAGE_SELF), time.monotonic()
            row["candidate"] = prepare_partition(folder / "graph/orthogroups_profiles_refined.txt",
                folder / "candidate", names, species, (q, t, scores), True,
                _expand_phylogeny_candidates, load_membership_constraints)
            after = resource.getrusage(resource.RUSAGE_SELF)
            row["candidate_resources"] = {"wall_s": time.monotonic() - candidate_started,
                "cpu_s": after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
                "peak_parent_process_rss_GiB": after.ru_maxrss / 1024 ** 2,
                "memory_scope": "Parent lifetime high-water mark, not an aggregate process-tree peak"}
            candidate = Path(row["candidate"]["candidate_partition"]["path"])
            constraints = Path(row["candidate"]["membership_constraints"]["path"])
            measurement_only = selective_experiment and selective_mode == "hmmer_score_only"
            exact_control = mode == arms[0] or measurement_only
            reference = case
            if recovery and measurement_only:
                reference = {
                    "candidate": projection_control["candidate"]["candidate_partition"],
                    "constraints": projection_control["candidate"]["membership_constraints"],
                    "golden": {Path(pin["path"]).name: pin for pin in projection_control["outputs"]},
                }
                reference["golden"]["species_tree.rooted.nwk"] = projection_control["rooted_species_tree"]
            if exact_control:
                row["candidate_control_matches"] = candidate_control_matches(
                    candidate, constraints, reference["candidate"], reference["constraints"])
            row["native_execution"] = timed_stage(native_command, folder, "native", source, env)
            if row["native_execution"]["returncode"]:
                raise RuntimeError("Native phylogeny replay failed")
            phylo = folder / "native/orthohmm_phylogeny"
            prediction = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
            row["coverage"] = validate_partition(candidate, prediction)
            row["outputs"] = [record(phylo / name) for name in PRODUCTS]
            if recovery:
                row["rooted_species_tree"] = record(phylo / "species_tree.rooted.nwk")
            native = replay_native(row, {"command": native_command})
            if (native["genes"], native["species"]) != UNIVERSE:
                raise ValueError("Wrong native universe")
            row["native_summary"] = native["summary"]
            if selective_experiment and mode != arms[0]:
                row["native_pair_turnover"] = native_pair_turnover(
                    Path(case["baseline_output"]) / "orthohmm_phylogeny/orthohmm_pairwise_orthologs_confidence.tsv",
                    phylo / "orthohmm_pairwise_orthologs_confidence.tsv")
            if exact_control:
                row["golden_matches"] = compare_golden(phylo, reference["golden"])
            row.update(score_orthobench(repo, prediction))
            actual, details, score_files = readback_scores(repo, "OrthoBench", row, native, folder, "weighting")
            if actual != row["scores"]:
                raise ValueError("Independent score readback differs")
            row.update(details=details, scoring_readback=score_files, metric_vector=endpoint_vector("OrthoBench", actual))
            for pin in (*checked, *native["checked"], *score_files):
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Runtime changed during evaluation")
            if mode == arms[0]:
                if normalized(row) != case["expected"]:
                    raise ValueError("Replayed control scores differ from admitted canonical control")
                control = normalized(row)
            elif measurement_only:
                if normalized(row) != (normalized(projection_control) if recovery else control):
                    raise ValueError("HMMER measurement-only scores differ from off")
                row["measurement_control_matches"] = True
                measurement_control = True
            else:
                row["advancement"] = advancement(control, normalized(row))
                if recovery and mode == "projection_off":
                    projection_control = row
                if selective_experiment and selective_mode != "off":
                    row["accuracy_interpretation"] = row["activation"]["applicability"]
                    if row["activation"]["completely_measured_queries"] == 0:
                        row["advancement"]["eligible"] = False
                        row["advancement"]["reason"] = "No completely measured queries; applicability fallback remains required"
                if recovery and selective_mode != "off":
                    row["matched_projection_advancement"] = advancement(normalized(projection_control), normalized(row))
            row["status"] = "independent_content_readback_verified"
        except Exception as error:
            row.update(status="failed", error_type=type(error).__name__, error=str(error))
        finally:
            row["total_incremental_wall_s"] = time.monotonic() - started
            write_json(output / "results.json", report)
        print(mode, row["status"], flush=True)
    report["status"] = ("finished_all_arms_verified" if all(
        row["status"] == "independent_content_readback_verified" for row in report["arms"])
        else "finished_with_explicit_gaps")
    report["host_after"] = host_activity()
    write_json(output / "results.json", report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands_parser = parser.add_subparsers(dest="command", required=True)
    preparation = commands_parser.add_parser("prepare")
    preparation.add_argument("--experiment", choices=("weighting", "alignment", "selective_hmmer", "projection_recovery"), default="weighting")
    for name in ("repo", "source", "output"):
        preparation.add_argument("--" + name, required=True, type=Path)
    execution = commands_parser.add_parser("run")
    execution.add_argument("--execution-plan", required=True, type=Path)
    args = parser.parse_args(argv)
    result = (prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve(), args.experiment)
              if args.command == "prepare" else run(args.execution_plan.resolve()))
    print(result["status"])
    return 0 if result["status"] in ("prepared_unrun", "finished_all_arms_verified") else 1


if __name__ == "__main__":
    raise SystemExit(main())

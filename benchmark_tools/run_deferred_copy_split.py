"""Regenerate native candidates and screen bounded copy deferral with a fixed tree."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from Bio import SeqIO

from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_orthobench_factorial import prepare_partition
from benchmark_tools.probe_leiden_boundary import saved_fingerprint
from benchmark_tools.report_hillclimb import readback_scores, replay_native
from benchmark_tools.run_held_tree_hillclimb import (
    require_control_scores, require_reuse, source_records, supplied_command, tree_contrast,
)
from benchmark_tools.run_membership_hillclimb import (
    PRODUCTS, argument, compare_golden, score_qfo, validate_partition, write_json,
)
from benchmark_tools.run_search_recovery import host_activity, runtime_records
from benchmark_tools.run_three_kingdoms_membership_supplement import (
    candidate_control_matches, timed_stage, validate_ownership,
)

PROTOCOL = "hillclimb_deferred_copy_split_20261013.json"
LABELS = ("control", "deferred500")


def regenerate_refinement(clusters, edges, species, names, preset, seed, hits=None):
    from orthohmm.refinement import refine_cluster_indices, resolve_refinement_profile
    empty = np.array([], dtype=np.int32)
    if hits is None:
        hits = (empty, empty, np.array([], dtype=np.float64))
    groups = refine_cluster_indices(clusters, *hits,
        edges[0], edges[1], species, rbnh_scores=edges[2], **resolve_refinement_profile(preset))
    with seed.open("x") as stream:
        for group in groups:
            stream.write(" ".join(sorted(names[i] for i in group)) + "\n")
    read_partition(seed, set(names))
    return groups


def refinement_changes(before, control, variant, species, maximum):
    def assignments(groups):
        labels = np.full(len(species), -1, dtype=np.int32)
        for i, group in enumerate(groups):
            if not group or np.any(labels[group] >= 0):
                raise ValueError("Duplicate or empty refined group")
            labels[group] = i
        if np.any(labels < 0):
            raise ValueError("Refinement lost input genes")
        return labels
    old, new = assignments(control), assignments(variant)
    assignments(before)
    counts = Counter()
    sizes = Counter()
    for group in before:
        copies = Counter(int(species[g]) for g in group)
        if len(set(old[group])) > 1 and len(set(new[group])) == 1:
            if len(group) > maximum or max(copies.values()) < 10:
                raise ValueError("Deferral exceeded its size/copy safeguard")
            counts["retained_pre_refinement_families"] += 1
            counts["genes_in_retained_families"] += len(group)
            sizes[str(len(group))] += 1
    return {**counts, "retained_family_size_counts": dict(sizes),
            "control_family_size_counts": dict(Counter(str(len(g)) for g in control)),
            "variant_family_size_counts": dict(Counter(str(len(g)) for g in variant))}


def prepare(repo, source, output, graph_plan, comparison, cpu):
    if output.exists() or isinstance(cpu, bool) or not 1 <= cpu <= 16:
        raise ValueError("Require fresh output and 1-16 CPUs")
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools"], check=True)
    graphs = json.loads(graph_plan.read_text())
    graph = graphs["graphs"]["profile_expanded"]
    partition = graphs["partitions"]["profile_expanded"]["file"]
    golden_refinement = graphs["partitions"]["profiles_refined"]["file"]
    factorial_path = repo / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
    factorial = json.loads(factorial_path.read_text())
    if graphs["partitions"]["candidate"]["file"] != factorial["candidate_arms"]["p1_c1"]["candidate_partition"]:
        raise ValueError("Native graph and candidate ancestry differ")
    numeric = factorial["numeric_checkpoint"]["manifest"]
    checkpoint = Path(numeric["path"]).parent
    report = json.loads(comparison.read_text())
    method = next(m for m in report["methods"] if m["key"] == "membership:control_direct")
    cells = [c for c in method["cells"] if c["dataset"] == "QfO"
             and c["status"] == "independent_content_readback_verified"]
    if len(cells) != 1:
        raise ValueError("Require one verified corrected QfO native control")
    cell = cells[0]
    check(cell["source_snapshot"])
    producer = json.loads(Path(cell["source_snapshot"]["path"]).read_text())
    check(producer["execution_plan"])
    previous = json.loads(Path(producer["execution_plan"]["path"]).read_text())
    original = next(c["command"] for c in previous["cases"]
                    if c["dataset"] == "QfO" and c["label"] == "control_direct")
    checkpoint_root = Path(argument(original, "--output-directory")) / "orthohmm_phylogeny"
    tree = checkpoint_root / "species_tree.rooted.nwk"
    payload = Path(graph["payload"])
    inputs = [record(graph_plan), record(comparison), record(factorial_path),
              partition, golden_refinement, cell["source_snapshot"], producer["execution_plan"],
              *cell["checked"], record(tree), *[record(p) for p in sorted(checkpoint.iterdir())],
              *[record(payload / name) for name in ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")]]
    for directory in ("alignments", "checkpoints", "gene_trees"):
        files = [p for p in sorted((checkpoint_root / directory).rglob("*")) if p.is_file()]
        if not files:
            raise ValueError("Missing native tree checkpoint directory")
        inputs.extend(record(p) for p in files)
    baseline = factorial["candidate_arms"]["p1_c1"]
    inputs.extend([baseline["candidate_partition"], baseline["membership_constraints"],
                   record(source / "benchmark_tools" / PROTOCOL)])
    inputs = list({item["path"]: item for item in inputs}.values())
    for item in inputs:
        check(item)
    nonce = hashlib.sha256(str(output).encode()).hexdigest()[:16]
    cases = []
    for i, label in enumerate(LABELS):
        folder = output / label
        command = supplied_command(original, source, folder, tree, cpu)
        working = folder / "candidate/orthohmm_working_res"
        for flag, value in (("--candidate-clusters", working / "orthohmm_edges_clustered.txt"),
                            ("--membership-constraints", working / "phylogeny_candidate_merges.json")):
            argument(command, flag)
            command[command.index(flag) + 1] = str(value)
        cases.append({"label": label, "output": str(folder), "command": command,
                      "scoring_id": f"deferred_{nonce}_{i}"})
    plan = {"status": "prepared_unrun", "repo": str(repo), "source_root": str(source),
            "output": str(output), "cpu": cpu, "inputs": inputs,
            "sources": [*source_records(source), record(source / "benchmark_tools/run_deferred_copy_split.slurm")],
            "runtime": runtime_records(), "tool_path": os.environ["PATH"],
            "source_revision": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
            "graph": graph, "partition": partition, "golden_refinement": golden_refinement,
            "checkpoint": str(checkpoint), "checkpoint_manifest": numeric,
            "fasta_directory": argument(original, "--fasta-directory"), "tree": record(tree),
            "baseline_candidate": baseline["candidate_partition"],
            "baseline_constraints": baseline["membership_constraints"],
            "baseline_candidate_parameters": baseline["expansion"]["parameters"],
            "golden": cell["products"], "cases": cases,
            "expected_control_scores": {k.split(": ", 1)[1]: v for k, v in cell["scores"].items()},
            "scope": "Development-exposed final-refinement/candidate change with fixed species tree; cached HMM graph and compatible raw gene trees, not fresh pipeline timing."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", plan)
    return plan


def construct(plan, report, destination):
    from orthohmm.accuracy import load_accuracy_checkpoint, read_index_clusters
    from orthohmm.orthohmm import _expand_phylogeny_candidates
    from benchmark_tools.replay_phylogeny import load_membership_constraints

    checkpoint = Path(plan["checkpoint"])
    report["checkpoint_audit"] = audit(checkpoint, plan["checkpoint_manifest"]["sha256"])
    names, species, queries, targets, scores = load_accuracy_checkpoint(checkpoint, verify=False)
    owners = {}
    for fasta in sorted(Path(plan["fasta_directory"]).glob("*.fasta")):
        for sequence in SeqIO.parse(fasta, "fasta"):
            if sequence.id in owners or not sequence.seq:
                raise ValueError("Duplicate or empty FASTA gene")
            owners[sequence.id] = fasta.name
    report["ownership"] = validate_ownership(names, species, owners)
    if report["ownership"]["genes"] != 984137 or report["ownership"]["species"] != 78:
        raise ValueError("Wrong corrected QfO input universe")
    payload = Path(plan["graph"]["payload"])
    if (payload / "gene_names.txt").read_text().splitlines() != names:
        raise ValueError("Native graph/checkpoint gene order differs")
    if saved_fingerprint(payload) != plan["graph"]["fingerprint"]:
        raise ValueError("Native graph fingerprint differs")
    clusters = read_index_clusters(plan["partition"]["path"], names)
    if len(clusters) != 306051:
        raise ValueError("Native profile-expanded partition count differs")
    read_partition(Path(plan["partition"]["path"]), set(names))
    edges = [np.load(payload / (name + ".npy"), mmap_mode="r", allow_pickle=False)
             for name in ("sources", "targets", "weights")]
    refined = {}
    for label, preset in zip(LABELS, ("default", "phylogeny_deferred")):
        started = time.monotonic()
        folder = destination / label
        folder.mkdir(exist_ok=False)
        seed = folder / "refined_seed.txt"
        groups = regenerate_refinement(clusters, edges, species, names, preset, seed)
        if label == "control" and record(seed)["sha256"] != plan["golden_refinement"]["sha256"]:
            raise ValueError("Default final refinement does not reproduce admitted native bytes")
        refined[label] = groups
        row = {"label": label, "status": "candidate_preparation_running", "scores": {}, "seed": record(seed)}
        report["arms"].append(row)
        write_json(destination / "results.json", report)
        row["candidate"] = prepare_partition(seed, folder / "candidate", names, species,
            (queries, targets, scores), True, _expand_phylogeny_candidates, load_membership_constraints)
        if row["candidate"]["expansion"]["parameters"] != plan["baseline_candidate_parameters"]:
            raise ValueError("Satellite candidate parameters changed")
        if label == "control":
            row["candidate_control_matches"] = candidate_control_matches(
                Path(row["candidate"]["candidate_partition"]["path"]),
                Path(row["candidate"]["membership_constraints"]["path"]),
                plan["baseline_candidate"], plan["baseline_constraints"])
        row.update(status="candidates_verified_unscored", candidate_wall_s=time.monotonic() - started)
        write_json(destination / "results.json", report)
    report["refinement_changes"] = refinement_changes(clusters, refined["control"], refined["deferred500"], species, 500)
    write_json(destination / "results.json", report)


def run(plan_path):
    identity = record(plan_path)
    plan = json.loads(plan_path.read_text())
    source, repo, output = (Path(plan[key]) for key in ("source_root", "repo", "output"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_deferred_copy_split.py":
        raise ValueError("Execute the frozen source")
    if tuple(c["label"] for c in plan["cases"]) != LABELS or (output / "results.json").exists():
        raise ValueError("Unexpected case inventory or existing attempt")
    evaluation_mode = plan.get("evaluation_mode", "fixed_species_tree")
    if evaluation_mode not in ("fixed_species_tree", "inferred_parent"):
        raise ValueError("Unknown deferred copy-split evaluation mode")
    inferred_parent = evaluation_mode == "inferred_parent"
    if inferred_parent and any(argument(c["command"], "--species-tree-mode") != "infer"
                               or "--species-tree" in c["command"] for c in plan["cases"]):
        raise ValueError("Own-tree follow-up requires infer mode without a supplied tree")
    checked = [identity, *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for item in checked:
        check(item)
    if runtime_records() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Changed runtime or CUDA visibility")
    report = {"status": "running", "execution_plan": identity, "source_revision": plan["source_revision"],
              "arms": [], "job_id": os.environ.get("SLURM_JOB_ID"), "scope": plan["scope"],
              "host_activity_before": host_activity(), "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", report)
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "PYTHONHASHSEED": "0",
           "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    try:
        if inferred_parent:
            from benchmark_tools.deferred_qfo_inferred import resolve
            resolved = resolve(plan, output)
            report["arms"] = resolved["arms"]
            report["parent_admission"] = resolved["receipt"]
            report["parent_readbacks"] = resolved["readbacks"]
            checked.extend(resolved["checked"])
            for case in plan["cases"]:
                Path(case["output"]).mkdir(exist_ok=False)
        else:
            construct(plan, report, output)
        for item in checked:
            check(item)
        control_valid = False
        for row, case in zip(report["arms"], plan["cases"]):
            if case["label"] != "control" and not control_valid:
                row["status"] = "not_run_control_gate_failed"
                continue
            started = time.monotonic()
            folder = Path(case["output"])
            try:
                row["status"] = "native_running"
                write_json(output / "results.json", report)
                row["native_execution"] = timed_stage(case["command"], folder, "native", source, env)
                if row["native_execution"]["returncode"]:
                    raise RuntimeError("Native phylogenetic replay failed")
                phylo = folder / "output/orthohmm_phylogeny"
                prediction = folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"
                row["coverage"] = validate_partition(Path(argument(case["command"], "--candidate-clusters")), prediction)
                row["outputs"] = [record(phylo / name) for name in PRODUCTS]
                native = replay_native(row, case)
                if native["genes"] != 984137 or native["species"] != 78:
                    raise ValueError("Wrong native input universe")
                if inferred_parent:
                    from benchmark_tools.deferred_qfo_inferred import require_inferred_summary
                    require_inferred_summary(native["summary"])
                elif case["label"] == "control":
                    require_reuse(native["summary"])
                if case["label"] == "control":
                    row["golden_matches"] = compare_golden(phylo, plan["golden"])
                row["tree_handoff"] = tree_contrast(plan["tree"]["path"], phylo / "species_tree.rooted.nwk")
                if not inferred_parent and row["tree_handoff"]["rooted_clade_distance"]:
                    raise ValueError("Fixed species-tree handoff changed")
                row["native_summary"] = native["summary"]
                row["status"] = "native_outputs_verified_scoring_pending"
                write_json(output / "results.json", report)
                row.update(score_qfo(repo, phylo, folder, case["scoring_id"]))
                scores, details, files = readback_scores(repo, "QfO", row, native, folder,
                                                        "ohmm_hillclimb_r1_" + case["scoring_id"])
                if case["label"] == "control":
                    require_control_scores(scores, plan["expected_control_scores"])
                for item in [*checked, *native["checked"], *files]:
                    check(item)
                row.update(scores=scores, details=details, scoring_readback=files,
                           status="independent_content_readback_verified", inputs_reverified_unchanged=True)
                if case["label"] == "control":
                    control_valid = True
            except Exception as error:
                row.update(status="failed", error_type=type(error).__name__, error=str(error))
            finally:
                row["native_and_scoring_wall_s"] = time.monotonic() - started
                write_json(output / "results.json", report)
        report["status"] = ("finished_all_arms_verified" if all(row["status"] == "independent_content_readback_verified"
                             for row in report["arms"]) else "finished_with_explicit_gaps")
    except Exception as error:
        report.update(status="failed_preparation_or_input_gate", error_type=type(error).__name__, error=str(error))
    finally:
        try:
            for item in checked:
                check(item)
        except Exception as error:
            report.update(status="failed_input_reverification", error_type=type(error).__name__, error=str(error))
        report["host_activity_after"] = host_activity()
        write_json(output / "results.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="mode", required=True)
    preparation = commands.add_parser("prepare")
    for name in ("repo", "source-root", "output", "graph-plan", "comparison"):
        preparation.add_argument("--" + name, type=Path, required=True)
    preparation.add_argument("--cpu", type=int, default=16)
    inferred = commands.add_parser("prepare-inferred")
    for name in ("repo", "source-root", "output", "parent-plan", "parent-results"):
        inferred.add_argument("--" + name, type=Path, required=True)
    inferred.add_argument("--cpu", type=int, default=16)
    inferred.add_argument("--parent-job", required=True)
    execution = commands.add_parser("run")
    execution.add_argument("--plan", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "prepare":
        result = prepare(args.repo.resolve(), args.source_root.resolve(), args.output.resolve(), args.graph_plan.resolve(),
                         args.comparison.resolve(), args.cpu)
    elif args.mode == "prepare-inferred":
        from benchmark_tools.deferred_qfo_inferred import prepare as prepare_inferred
        result = prepare_inferred(args.repo.resolve(), args.source_root.resolve(), args.output.resolve(), args.parent_plan.resolve(),
                                  args.parent_results.resolve(), args.cpu, args.parent_job)
    else:
        result = run(args.plan.resolve())
    print(json.dumps({"status": result["status"]}))
    return 0 if result["status"] in ("prepared_unrun", "prepared_waiting_for_terminal_parent", "finished_all_arms_verified") else 1


if __name__ == "__main__":
    raise SystemExit(main())

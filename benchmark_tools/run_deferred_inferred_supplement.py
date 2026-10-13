"""Evaluate bounded copy deferral with candidate-specific inferred species trees."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from Bio import SeqIO
from orthohmm.files import fetch_fasta_files

from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_orthobench_factorial import prepare_partition
from benchmark_tools.replay_high_sensitivity import production_refinement_hits
from benchmark_tools.report_hillclimb import readback_scores, replay_native
from benchmark_tools.run_deferred_copy_split import regenerate_refinement, refinement_changes
from benchmark_tools.run_held_tree_hillclimb import source_records
from benchmark_tools.run_membership_hillclimb import (
    PRODUCTS, argument, compare_golden, derive_command, original_commands, scoring_records,
    score_orthobench, score_three_kingdoms, validate_partition, write_json,
)
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_three_kingdoms_membership_supplement import (
    candidate_control_matches, timed_stage, validate_ownership,
)

DATASETS = ("OrthoBench", "ThreeKingdoms")
LABELS = ("control", "deferred500")
PROTOCOL = "hillclimb_deferred_inferred_supplement_20261013.json"


def load_native_edges(path, names, count):
    from orthohmm.helpers import IndexedEdges
    if isinstance(count, bool) or not isinstance(count, int) or count <= 0 or len(set(names)) != len(names):
        raise ValueError("Require positive native count and unique gene table")
    index = {name: i for i, name in enumerate(names)}
    sources, targets = np.empty(count, dtype=np.int32), np.empty(count, dtype=np.int32)
    weights = np.empty(count, dtype=np.float64)
    seen, previous = 0, -1
    with path.open() as stream:
        for row in csv.reader(stream, delimiter="\t"):
            if len(row) != 3 or seen == count:
                raise ValueError("Native edge row/count differs")
            try:
                left, right, weight = index[row[0]], index[row[1]], float(row[2])
            except (KeyError, ValueError) as error:
                raise ValueError("Unknown native gene or malformed weight") from error
            packed = left * len(names) + right
            if left >= right or packed <= previous or not math.isfinite(weight) or weight <= 0:
                raise ValueError("Noncanonical, duplicate, unordered or invalid native edge")
            sources[seen], targets[seen], weights[seen] = left, right, weight
            seen += 1
            previous = packed
    if seen != count:
        raise ValueError("Native edge count differs")
    return IndexedEdges(names, sources, targets, weights)


def selected_datasets(names):
    names = tuple(names)
    if not names or tuple(d for d in DATASETS if d in names) != names:
        raise ValueError("Require a nonempty ordered unique subset of the supplement datasets")
    return names


def native_control_paths(dataset, native_output, original):
    if dataset == "OrthoBench":
        working = native_output / "orthohmm_working_res"
        return (native_output / "orthohmm_phylogeny", working / "phylogeny_candidate_superfamilies.txt",
                working / "phylogeny_candidate_merges.json", native_output)
    if dataset != "ThreeKingdoms":
        raise ValueError("Unknown inferred supplement dataset")
    baseline = Path(argument(original, "--output-directory"))
    return (baseline / "orthohmm_phylogeny", Path(argument(original, "--candidate-clusters")),
            Path(argument(original, "--membership-constraints")), baseline)


def prepare(repo, source, output, comparison, cpu, datasets=DATASETS):
    if output.exists() or isinstance(cpu, bool) or not 1 <= cpu <= 16:
        raise ValueError("Require fresh namespace and 1-16 CPUs")
    datasets = selected_datasets(datasets)
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools"], check=True)
    comparison_report = json.loads(comparison.read_text())
    method = next(m for m in comparison_report["methods"] if m["key"] == "membership:control_direct")
    cases, inputs = [], [record(comparison), record(source / "benchmark_tools" / PROTOCOL)]
    for dataset in datasets:
        cell = next(c for c in method["cells"] if c["dataset"] == dataset)
        if cell["status"] != "independent_content_readback_verified":
            raise ValueError("Require independently validated control cell")
        original, _, evidence = original_commands(repo, dataset)
        baseline = Path(argument(original, "--output-directory"))
        if dataset == "OrthoBench":
            readback_path = repo / "benchmark_tools/results/hillclimb_controls_readback_v1_20261012/report.json"
            readback = json.loads(readback_path.read_text())
            if readback["fresh_orthobench_control"]["status"] != "independent_content_readback_verified":
                raise ValueError("Fresh OrthoBench control has no independent readback")
            capture = readback["captures"]["full_orthobench"]["snapshot"]
            check(capture)
            producer_path = Path(capture["path"])
            producer = json.loads(producer_path.read_text())
            control = next(a for a in producer["arms"] if a["label"] == "control100")
            if control["status"] != "scored_development_result":
                raise ValueError("Fresh cap100 OrthoBench control is not scored")
            native_output = Path(control["native"]["metrics"]["metadata"]["output_directory"])
            metrics_path = native_output.parent / "metrics.json"
            inputs.extend([record(readback_path), capture, producer["execution_plan"],
                           *readback["fresh_orthobench_control"]["checked"]])
        else:
            native_output = baseline
            metrics_path = repo / "three_kingdoms/results/parity_20260907/orthohmm_phylogeny_satellite_v2/metrics.json"
        metrics = json.loads(metrics_path.read_text())
        metadata = metrics["metadata"]
        expected = {"accuracy_profile": "high_sensitivity", "clustering": "leiden", "cpm_resolution": .1,
                    "substitution_matrix": "BLOSUM62", "evalue_threshold": 1e-4, "phylogeny_candidates": "satellite_v2",
                    "search_mode": "builtin", "search_max_candidates_per_query": 100, "leiden_seed": 4}
        if metrics["status"] != "complete" or any(metadata.get(k) != v for k, v in expected.items()):
            raise ValueError("Wrong native HMM/profile graph settings")
        fasta = Path(argument(original, "--fasta-directory"))
        if Path(metadata["fasta_directory"]).resolve() != fasta.resolve():
            raise ValueError("Native graph and control FASTA directories differ")
        working = native_output / "orthohmm_working_res"
        checkpoint = working / "high_sensitivity_checkpoint"
        if (Path(metadata["output_directory"]).resolve() != native_output.resolve()
                or Path(metadata["high_sensitivity_checkpoint"]).resolve() != checkpoint.resolve()):
            raise ValueError("Native graph and checkpoint provenance directories differ")
        tools = [record(shutil.which(argument(original, flag)) or argument(original, flag))
                 for flag in ("--aligner", "--tree-builder")]
        phylo, control_candidate, control_constraints, control_checkpoint = native_control_paths(dataset, native_output, original)
        golden = {name: record(phylo / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")}
        baseline_candidate = record(control_candidate)
        baseline_constraints = record(control_constraints)
        dataset_inputs = [*evidence, record(metrics_path), *cell["checked"], cell["source_snapshot"],
            *tools, *golden.values(), baseline_candidate, baseline_constraints,
            record(working / "orthohmm_edges.txt"), *[record(p) for p in sorted(checkpoint.iterdir())],
            *[record(fasta / name) for name in sorted(fetch_fasta_files(str(fasta)))], *scoring_records(repo, dataset)]
        for directory in ("alignments", "checkpoints", "gene_trees"):
            files = [p for p in sorted((phylo / directory).rglob("*")) if p.is_file()]
            if not files:
                raise ValueError("Missing baseline native tree cache")
            dataset_inputs.extend(record(p) for p in files)
        inputs.extend(dataset_inputs)
        arms = []
        for label in LABELS:
            folder = output / dataset / label
            command = derive_command(original, source, folder, "high_confidence_pair", cpu)
            command[command.index("--checkpoint-source") + 1] = str(control_checkpoint)
            for flag, value in (("--candidate-clusters", folder / "candidate/orthohmm_working_res/orthohmm_edges_clustered.txt"),
                                ("--membership-constraints", folder / "candidate/orthohmm_working_res/phylogeny_candidate_merges.json")):
                command[command.index(flag) + 1] = str(value)
            arms.append({"label": label, "command": command, "output": str(folder)})
        cases.append({"dataset": dataset, "fasta_directory": str(fasta), "checkpoint": str(checkpoint),
            "checkpoint_manifest": record(checkpoint / "manifest.json"),
            "native_graph": record(working / "orthohmm_edges.txt"), "native_edges": metrics["counts"]["network_edges"],
            "expected_genes": {"OrthoBench": 251378, "ThreeKingdoms": 443217}[dataset], "expected_species": 12,
            "baseline_candidate": baseline_candidate, "baseline_constraints": baseline_constraints,
            "golden": golden, "expected_control_scores": (
                readback["fresh_orthobench_control"]["scores"] if dataset == "OrthoBench" else cell["scores"]), "arms": arms})
    inputs = list({item["path"]: item for item in inputs}.values())
    for item in inputs:
        check(item)
    result = {"status": "prepared_unrun", "repo": str(repo), "source_root": str(source), "output": str(output),
        "source_revision": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "cpu": cpu, "datasets": list(datasets), "inputs": inputs, "sources": [*source_records(source),
            record(source / "benchmark_tools/run_deferred_inferred_supplement.slurm")],
        "runtime": runtime(), "tool_path": os.environ["PATH"], "cases": cases,
        "scope": "Development-exposed lossless native HMM/profile graph reconstruction with candidate-specific inferred phylogeny; not fresh full-pipeline timing."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def construct(case, report, folder):
    from orthohmm.accuracy import load_accuracy_checkpoint, read_index_clusters
    from orthohmm.externals import execute_leiden
    from orthohmm.orthohmm import _expand_phylogeny_candidates
    from benchmark_tools.replay_phylogeny import load_membership_constraints

    checkpoint = Path(case["checkpoint"])
    report["checkpoint_audit"] = audit(checkpoint, case["checkpoint_manifest"]["sha256"])
    names, species, queries, targets, scores = load_accuracy_checkpoint(checkpoint, verify=False)
    owners = {}
    fasta_root = Path(case["fasta_directory"])
    for name in sorted(fetch_fasta_files(str(fasta_root))):
        fasta = fasta_root / name
        for sequence in SeqIO.parse(fasta, "fasta"):
            if sequence.id in owners or not sequence.seq:
                raise ValueError("Duplicate or empty FASTA gene")
            owners[sequence.id] = fasta.name
    report["ownership"] = validate_ownership(names, species, owners)
    if report["ownership"]["genes"] != case["expected_genes"] or report["ownership"]["species"] != case["expected_species"]:
        raise ValueError("Wrong dataset gene/species universe")
    edges = load_native_edges(Path(case["native_graph"]["path"]), names, case["native_edges"])
    graph = folder / "graph"
    (graph / "orthohmm_working_res").mkdir(parents=True)
    execute_leiden(.1, str(graph), edges=edges, include_isolates=True, seed=4)
    partition = graph / "orthohmm_working_res/orthohmm_edges_clustered.txt"
    read_partition(partition, set(names))
    clusters = read_index_clusters(str(partition), names)
    report["graph_partition"] = record(partition)
    hits = production_refinement_hits(queries, targets, scores, species)
    report["refinement_directed_hits"] = len(hits[2])
    refined = {}
    for arm, preset in zip(case["arms"], ("default", "phylogeny_deferred")):
        arm_folder = Path(arm["output"])
        arm_folder.mkdir()
        seed = arm_folder / "refined_seed.txt"
        groups = regenerate_refinement(clusters, (edges.sources, edges.targets, edges.weights), species, names,
                                       preset, seed, hits=hits)
        refined[arm["label"]] = groups
        candidate = prepare_partition(seed, arm_folder / "candidate", names, species, (queries, targets, scores),
                                      True, _expand_phylogeny_candidates, load_membership_constraints)
        row = {"label": arm["label"], "status": "candidates_verified_unscored", "scores": {}, "candidate": candidate}
        if arm["label"] == "control":
            row["candidate_control_matches"] = candidate_control_matches(
                Path(candidate["candidate_partition"]["path"]), Path(candidate["membership_constraints"]["path"]),
                case["baseline_candidate"], case["baseline_constraints"])
        report["arms"].append(row)
        write_json(folder / "results.json", report)
    report["refinement_changes"] = refinement_changes(clusters, refined["control"], refined["deferred500"], species, 500)
    write_json(folder / "results.json", report)


def run(path):
    identity = record(path)
    plan = json.loads(path.read_text())
    source, repo, output = (Path(plan[k]) for k in ("source_root", "repo", "output"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_deferred_inferred_supplement.py":
        raise ValueError("Execute frozen source")
    datasets = selected_datasets(plan.get("datasets", DATASETS))
    if tuple(c["dataset"] for c in plan["cases"]) != datasets or (output / "results.json").exists():
        raise ValueError("Unexpected dataset inventory or existing attempt")
    for case in plan["cases"]:
        if tuple(a["label"] for a in case["arms"]) != LABELS or any(
            argument(a["command"], "--species-tree-mode") != "infer" or "--species-tree" in a["command"]
            for a in case["arms"]):
            raise ValueError("Require both matched arms and candidate-specific inferred trees")
    checked = [identity, *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for item in checked:
        check(item)
    if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Changed runtime or CUDA visibility")
    result = {"status": "running", "execution_plan": identity, "source_revision": plan["source_revision"],
              "datasets": [], "job_id": os.environ.get("SLURM_JOB_ID"), "scope": plan["scope"],
              "host_activity_before": host_activity(), "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", result)
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "PYTHONHASHSEED": "0",
           "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    for case in plan["cases"]:
        folder = output / case["dataset"]
        folder.mkdir()
        report = {"dataset": case["dataset"], "status": "running", "arms": []}
        result["datasets"].append(report)
        try:
            construct(case, report, folder)
            generated = [item for row in report["arms"] for item in
                         (row["candidate"]["candidate_partition"], row["candidate"]["membership_constraints"], row["candidate"]["seed_partition"])]
            control_valid = False
            for row, arm in zip(report["arms"], case["arms"]):
                if arm["label"] != "control" and not control_valid:
                    row["status"] = "not_run_control_gate_failed"
                    continue
                started = time.monotonic()
                try:
                    for item in [*checked, *generated]:
                        check(item)
                    row["status"] = "native_running"
                    write_json(folder / "results.json", report)
                    arm_folder = Path(arm["output"])
                    row["native_execution"] = timed_stage(arm["command"], arm_folder, "native", source, env)
                    if row["native_execution"]["returncode"]:
                        raise RuntimeError("Native inferred phylogeny failed")
                    phylo = arm_folder / "output/orthohmm_phylogeny"
                    prediction = arm_folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"
                    row["coverage"] = validate_partition(Path(argument(arm["command"], "--candidate-clusters")), prediction)
                    row["outputs"] = [record(phylo / name) for name in PRODUCTS]
                    native = replay_native(row, arm)
                    if native["genes"] != case["expected_genes"] or native["species"] != case["expected_species"]:
                        raise ValueError("Wrong native universe")
                    if arm["label"] == "control":
                        row["golden_matches"] = compare_golden(phylo, case["golden"])
                    row["native_summary"] = native["summary"]
                    row.update(score_orthobench(repo, prediction) if case["dataset"] == "OrthoBench"
                               else score_three_kingdoms(repo, prediction, arm_folder))
                    scores, details, files = readback_scores(repo, case["dataset"], row, native, arm_folder, "unused")
                    if arm["label"] == "control" and {
                        key.split(": ", 1)[1]: value for key, value in case["expected_control_scores"].items()} != scores:
                        raise ValueError("Deterministic native control scores differ")
                    for item in [*checked, *generated, *native["checked"], *files]:
                        check(item)
                    row.update(scores=scores, details=details, scoring_readback=files,
                               status="independent_content_readback_verified", inputs_reverified_unchanged=True)
                    if arm["label"] == "control":
                        control_valid = True
                except Exception as error:
                    row.update(status="failed", error_type=type(error).__name__, error=str(error))
                finally:
                    row["native_and_scoring_wall_s"] = time.monotonic() - started
                    write_json(folder / "results.json", report)
            report["status"] = ("finished_all_arms_verified" if all(r["status"] == "independent_content_readback_verified"
                                for r in report["arms"]) else "finished_with_explicit_gaps")
        except Exception as error:
            report.update(status="failed_preparation_or_gate", error_type=type(error).__name__, error=str(error))
        write_json(folder / "results.json", report)
        write_json(output / "results.json", result)
    result["status"] = ("finished_all_datasets_verified" if all(r["status"] == "finished_all_arms_verified"
                        for r in result["datasets"]) else "finished_with_explicit_gaps")
    try:
        for item in checked:
            check(item)
    except Exception as error:
        result.update(status="failed_input_reverification", error_type=type(error).__name__, error=str(error))
    result["host_activity_after"] = host_activity()
    write_json(output / "results.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="mode", required=True)
    prep = commands.add_parser("prepare")
    for name in ("repo", "source-root", "output", "comparison"):
        prep.add_argument("--" + name, type=Path, required=True)
    prep.add_argument("--cpu", type=int, default=16)
    prep.add_argument("--datasets", nargs="+", choices=DATASETS, default=list(DATASETS))
    execution = commands.add_parser("run")
    execution.add_argument("--plan", type=Path, required=True)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source_root.resolve(), args.output.resolve(), args.comparison.resolve(), args.cpu, args.datasets)
              if args.mode == "prepare" else run(args.plan.resolve()))
    print(json.dumps({"status": result["status"]}))
    return 0 if result["status"] in ("prepared_unrun", "finished_all_datasets_verified") else 1


if __name__ == "__main__":
    raise SystemExit(main())

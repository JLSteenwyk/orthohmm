"""Compare complete downstream inference for arrival/canonical hit block order."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from Bio import SeqIO
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.audit_hit_arrival_order import block_layout
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_orthobench_factorial import prepare_partition
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import (
    DATASETS, PRODUCTS, argument, compare_golden, original_commands,
    score_orthobench, score_qfo, score_three_kingdoms, scoring_records,
    validate_partition, write_json,
)
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import (
    candidate_control_matches, replay_command, sources, timed_stage, validate_ownership,
)
from benchmark_tools.validate_profile_runtime import require_profile_runtime
from orthohmm.accuracy import load_accuracy_checkpoint, write_accuracy_checkpoint
from orthohmm.files import fetch_fasta_files
from orthohmm.orthohmm import _expand_phylogeny_candidates

PLAN = "hillclimb_canonical_downstream_20261014.json"
LABELS = ("arrival", "canonical")
GENES = {"OrthoBench": 251378, "QfO": 984137, "ThreeKingdoms": 443217}


def filename_block_order(species, queries, targets, classes):
    starts, ends, keys = block_layout(species, queries, targets)
    base = int(np.max(species)) + 1
    labels = set(int(value) for value in np.unique(species))
    if set(classes) != labels or len(set(classes.values())) != len(classes):
        raise ValueError("Require bijective complete species/filename classes")
    order = sorted(range(len(keys)), key=lambda i: (
        classes[int(keys[i]) // base], classes[int(keys[i]) % base]))
    return starts, ends, order


def ordered_checkpoint(checkpoint, expected_sha, owners, output):
    admission = audit(checkpoint, expected_sha)
    names, species, queries, targets, scores = load_accuracy_checkpoint(checkpoint, verify=False)
    ownership = validate_ownership(names, species, owners)
    starts, ends, order = filename_block_order(species, queries, targets, ownership["classes"])
    arrays = tuple(np.concatenate([values[starts[i]:ends[i]] for i in order])
                   for values in (queries, targets, scores))
    # Compare each original block with its destination, preserving score bits and row order.
    offset = 0
    for i in order:
        count = int(ends[i] - starts[i])
        for old, new in zip((queries, targets, scores), arrays):
            if old[starts[i]:ends[i]].tobytes() != new[offset:offset + count].tobytes():
                raise ValueError("Canonical block changed an accepted hit")
        offset += count
    if offset != len(scores):
        raise ValueError("Canonical block ordering lost hits")
    result = write_accuracy_checkpoint(str(output), names, species, *arrays)
    manifest = record(result / "manifest.json")
    return {"checkpoint": str(result), "manifest": manifest,
            "audit": audit(result, manifest["sha256"]), "original_audit": admission,
            "ownership": ownership, "species_pair_blocks": len(order),
            "whole_block_multiset_and_internal_rows_preserved": True}


def baseline(repo, dataset):
    original, _, evidence = original_commands(repo, dataset)
    native = Path(argument(original, "--output-directory"))
    if dataset == "OrthoBench":
        path = repo / "benchmarks/results/full_search_hillclimb_r1_v3_20261012/OrthoBench/results.json"
        report = json.loads(path.read_text())
        row = next(row for row in report["arms"] if row["label"] == "control100")
        if report["status"] != "finished_all_arms_scored" or row["status"] != "scored_development_result":
            raise ValueError("Completed fresh OrthoBench control required")
        check(report["execution_plan"])
        plan = json.loads(Path(report["execution_plan"]["path"]).read_text())
        case = next(case for case in plan["cases"] if case["dataset"] == dataset and case["label"] == "control100")
        native = Path(case["output"]) / "output"
        checkpoint = native / "orthohmm_working_res/high_sensitivity_checkpoint"
        checkpoint_pin = row["native"]["numeric_checkpoint"]["manifest"]
        candidate = native / "orthohmm_working_res/phylogeny_candidate_superfamilies.txt"
        evidence.extend([record(path), report["execution_plan"]])
    elif dataset == "QfO":
        path = repo / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
        report = json.loads(path.read_text())
        checkpoint_pin = report["numeric_checkpoint"]["manifest"]
        checkpoint = Path(checkpoint_pin["path"]).parent
        candidate = Path(argument(original, "--candidate-clusters"))
        evidence.append(record(path))
    else:
        path = native.parent / "metrics.json"
        metadata = json.loads(path.read_text())["metadata"]
        checkpoint = Path(metadata["high_sensitivity_checkpoint"])
        checkpoint_pin = record(checkpoint / "manifest.json")
        candidate = Path(argument(original, "--candidate-clusters"))
    check(checkpoint_pin)
    return {"dataset": dataset, "fasta_directory": argument(original, "--fasta-directory"),
            "baseline_output": str(native), "checkpoint": str(checkpoint),
            "checkpoint_manifest": checkpoint_pin, "evidence": evidence,
            "candidate": record(candidate),
            "constraints": record(candidate.parent / "phylogeny_candidate_merges.json"),
            "golden": {name: record(native / "orthohmm_phylogeny" / name)
                       for name in (*PRODUCTS, "species_tree.rooted.nwk")},
            "aligner": argument(original, "--aligner"), "tree_builder": argument(original, "--tree-builder")}


def prepare(repo, source, output):
    if output.exists() or source == repo:
        raise ValueError("Require new namespace and frozen source checkout")
    plan_path = source / "benchmark_tools" / PLAN
    plan = json.loads(plan_path.read_text())
    if tuple(plan["arms"]) != LABELS or tuple(plan["datasets"]) != DATASETS:
        raise ValueError("Unexpected prespecified case inventory")
    cases, inputs = [], [record(plan_path)]
    for dataset in DATASETS:
        case = baseline(repo, dataset)
        fasta = Path(case["fasta_directory"])
        case["fastas"] = [record(fasta / name) for name in fetch_fasta_files(str(fasta))]
        if len(case["fastas"]) != (78 if dataset == "QfO" else 12):
            raise ValueError("Unexpected FASTA inventory")
        case["scorers"] = scoring_records(repo, dataset)
        cases.append(case)
        inputs.extend([*case["evidence"], case["candidate"], case["constraints"],
                       *case["golden"].values(), *case["fastas"], *case["scorers"],
                       *[record(p) for p in sorted(Path(case["checkpoint"]).iterdir()) if p.is_file()]])
        phylo = Path(case["baseline_output"]) / "orthohmm_phylogeny"
        for directory in ("alignments", "gene_trees", "checkpoints", "species_tree_inference"):
            inputs.extend(record(p) for p in sorted((phylo / directory).rglob("*")) if p.is_file())
        inputs.extend(record(shutil.which(case[name]) or case[name]) for name in ("aligner", "tree_builder"))
    inputs = list({pin["path"]: pin for pin in inputs}.values())
    source_pins = sources(source) + [record(source / "benchmark_tools/run_canonical_downstream.slurm")]
    scientific_runtime = runtime()
    for pin in (*inputs, *source_pins, *scientific_runtime["files"]):
        check(pin)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    for pin in source_pins:
        relative = Path(pin["path"]).relative_to(source)
        if relative.suffix == ".so":
            continue
        blob = subprocess.check_output(["git", "show", revision + ":" + str(relative)], cwd=repo)
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Frozen source differs from committed revision: " + str(relative))
    original_source = repo / "benchmarks/work/full_search_hillclimb_r1_v3_20261012"
    for name in ("pair_align.so", "hmm_viterbi.so", "kmer_prefilter.so"):
        relative = Path("orthohmm/search/csrc") / name
        if record(source / relative)["sha256"] != record(original_source / relative)["sha256"]:
            raise ValueError("Require original unchanged CPU libraries")
    result = {"status": "prepared_unrun", "repo": str(repo), "source_root": str(source),
              "source_revision": revision, "output": str(output), "cpu": 16,
              "inputs": inputs, "sources": source_pins, "runtime": scientific_runtime,
              "profile_runtime": require_profile_runtime(source), "tool_path": os.environ["PATH"],
              "cases": cases, "scope": plan["timing_scope"], "defaults_promoted": False}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def commands(source, plan, case, folder, checkpoint, sha):
    execution = {"checkpoint": str(checkpoint), "plan": {"checkpoint_manifest_sha256": sha},
                 "fasta_directory": case["fasta_directory"], "cpu": plan["cpu"]}
    graph = replay_command(source, execution, folder, .1)
    working = folder / "candidate/orthohmm_working_res"
    native = [sys.executable, str(source / "benchmark_tools/replay_phylogeny.py"),
              "--fasta-directory", case["fasta_directory"],
              "--candidate-clusters", str(working / "orthohmm_edges_clustered.txt"),
              "--membership-constraints", str(working / "phylogeny_candidate_merges.json"),
              "--output-directory", str(folder / "native"), "--json", str(folder / "native_metrics.json"),
              "--cpu", str(plan["cpu"]), "--species-tree-mode", "infer",
              "--aligner", case["aligner"], "--tree-builder", case["tree_builder"],
              "--root-rule", "species_overlap", "--pair-rule", "positive_paralogy",
              "--species-tree-rooting", "min_variance", "--gene-tree-support-mode", "legacy",
              "--satellite-membership-rule", "high_confidence_pair",
              "--checkpoint-source", case["baseline_output"]]
    return graph, native


def run(plan_path):
    plan = json.loads(plan_path.read_text())
    source, repo, output = (Path(plan[key]) for key in ("source_root", "repo", "output"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_canonical_downstream.py":
        raise ValueError("Execute the frozen runner")
    if (output / "results.json").exists():
        raise FileExistsError("Preserve every previous attempt")
    checked = [record(plan_path), *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in checked:
        check(pin)
    if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Changed runtime or CUDA visibility")
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "PYTHONHASHSEED": "0",
           "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    result = {"status": "running", "execution_plan": record(plan_path), "job_id": os.environ.get("SLURM_JOB_ID"),
              "source_revision": plan["source_revision"], "arms": [], "scope": plan["scope"],
              "host_activity_before": host_activity(), "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", result)
    nonce = hashlib.sha256(str(output).encode()).hexdigest()[:16]
    for case in plan["cases"]:
        dataset = case["dataset"]
        owners = {}
        for pin in case["fastas"]:
            for sequence in SeqIO.parse(pin["path"], "fasta"):
                if sequence.id in owners or not sequence.seq:
                    raise ValueError("Duplicate or empty FASTA gene")
                owners[sequence.id] = Path(pin["path"]).name
        if len(owners) != GENES[dataset]:
            raise ValueError("Unexpected full gene universe")
        control_valid = False
        for label in LABELS:
            row = {"dataset": dataset, "label": label, "status": "running", "scores": {}}
            result["arms"].append(row)
            write_json(output / "results.json", result)
            if label == "canonical" and not control_valid:
                row["status"] = "not_run_failed_control_gate"
                write_json(output / "results.json", result)
                continue
            folder = output / dataset / label
            folder.mkdir(parents=True, exist_ok=False)
            started = time.monotonic()
            try:
                for pin in checked:
                    check(pin)
                checkpoint, sha = Path(case["checkpoint"]), case["checkpoint_manifest"]["sha256"]
                if label == "canonical":
                    row["ordering"] = ordered_checkpoint(checkpoint, sha, owners, folder / "ordered")
                    checkpoint, sha = Path(row["ordering"]["checkpoint"]), row["ordering"]["manifest"]["sha256"]
                else:
                    row["original_checkpoint_audit"] = audit(checkpoint, sha)
                names, species, q, t, scores = load_accuracy_checkpoint(checkpoint, verify=True)
                row["ownership"] = validate_ownership(names, species, owners)
                graph, native_command = commands(source, plan, case, folder, checkpoint, sha)
                write_json(folder / "preflight.json", {"graph": graph, "native": native_command,
                           "checkpoint": str(checkpoint), "checkpoint_manifest_sha256": sha})
                row["graph_execution"] = timed_stage(graph, folder, "graph", source, env)
                if row["graph_execution"]["returncode"]:
                    raise RuntimeError("Graph/profile replay failed")
                graph_report = json.loads((folder / "graph_metrics.json").read_text())
                for stage in graph_report["stages"]:
                    check(stage["output"])
                    read_partition(Path(stage["output"]["path"]), set(names))
                before = resource.getrusage(resource.RUSAGE_SELF)
                candidate_started = time.monotonic()
                row["candidate"] = prepare_partition(folder / "graph/orthogroups_profiles_refined.txt",
                    folder / "candidate", names, species, (q, t, scores), True,
                    _expand_phylogeny_candidates, load_membership_constraints)
                after = resource.getrusage(resource.RUSAGE_SELF)
                row["candidate_resources"] = {"wall_s": time.monotonic() - candidate_started,
                    "cpu_s": after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
                    "peak_process_rss_GiB": after.ru_maxrss / 1024 ** 2}
                candidate = Path(row["candidate"]["candidate_partition"]["path"])
                constraints = Path(row["candidate"]["membership_constraints"]["path"])
                if label == "arrival":
                    row["candidate_control_matches"] = candidate_control_matches(
                        candidate, constraints, case["candidate"], case["constraints"])
                row["native_execution"] = timed_stage(native_command, folder, "native", source, env)
                if row["native_execution"]["returncode"]:
                    raise RuntimeError("Native phylogeny replay failed")
                phylo = folder / "native/orthohmm_phylogeny"
                prediction = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
                row["coverage"] = validate_partition(candidate, prediction)
                row["outputs"] = [record(phylo / name) for name in PRODUCTS]
                native = replay_native(row, {"command": native_command})
                if native["genes"] != GENES[dataset] or native["species"] != (78 if dataset == "QfO" else 12):
                    raise ValueError("Wrong native universe")
                row["native_summary"] = native["summary"]
                if label == "arrival":
                    row["golden_matches"] = compare_golden(phylo, case["golden"])
                    control_valid = True
                row["status"] = "native_verified_scoring_pending"
                write_json(output / "results.json", result)
                if dataset == "OrthoBench":
                    row.update(score_orthobench(repo, prediction))
                elif dataset == "ThreeKingdoms":
                    row.update(score_three_kingdoms(repo, prediction, folder))
                else:
                    row.update(score_qfo(repo, phylo, folder, f"canon_{nonce}_{label}"))
                actual, details, score_files = readback_scores(repo, dataset, row, native, folder,
                    f"ohmm_hillclimb_r1_canon_{nonce}_{label}")
                if actual != row["scores"]:
                    raise ValueError("Score readback differs")
                row.update(details=details, scoring_readback=score_files,
                           metric_vector=endpoint_vector(dataset, actual))
                for pin in (*checked, *native["checked"], *score_files):
                    check(pin)
                if runtime() != plan["runtime"]:
                    raise ValueError("Runtime changed during evaluation")
                row["status"] = "independent_content_readback_verified"
            except Exception as error:
                row.update(status="failed", error_type=type(error).__name__, error=str(error))
                control_valid = False
            finally:
                row["total_incremental_wall_s"] = time.monotonic() - started
                write_json(output / "results.json", result)
            print(dataset, label, row["status"], flush=True)
    result["status"] = ("finished_all_arms_verified" if all(
        row["status"] == "independent_content_readback_verified" for row in result["arms"])
        else "finished_with_explicit_gaps")
    result["host_activity_after"] = host_activity()
    write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    preparation = sub.add_parser("prepare")
    for name in ("repo", "source", "output"):
        preparation.add_argument("--" + name, type=Path, required=True)
    execution = sub.add_parser("run")
    execution.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        report = prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve())
    else:
        report = run(args.execution_plan.resolve())
    print(report["status"])
    raise SystemExit(0 if report["status"] in {"prepared_unrun", "finished_all_arms_verified"} else 1)

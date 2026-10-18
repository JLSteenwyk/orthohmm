"""Run the prespecified complete-proteome OrthoBench profile-weighting pilots."""

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
from orthohmm.accuracy import load_accuracy_checkpoint
from orthohmm.orthohmm import _expand_phylogeny_candidates


ARMS = ("none", "henikoff", "henikoff_cluster", "species_henikoff_cluster")
PROTOCOL = "hmm_evidence_campaign_20261018.json"
RETAINED = "benchmarks/results/canonical_downstream_r1_v1_20261014/results.json"
UNIVERSE = (251378, 12)


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


def prepare(repo, source, output):
    if output.exists() or source == repo:
        raise ValueError("Require fresh output and a frozen scientific checkout")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
    protocol = source / "benchmark_tools" / PROTOCOL
    spec = json.loads(protocol.read_text())
    if tuple(arm["label"] for arm in spec["first_experiment"]["arms"]) != ARMS:
        raise ValueError("Unexpected prespecified weighting arms")
    case = selected_baseline(repo)
    if len(case["fastas"]) != 12:
        raise ValueError("Require all12 OrthoBench proteomes")
    inputs = [record(protocol), *case["evidence"], *case["fastas"], case["candidate"],
              case["constraints"], *case["golden"].values(), *scoring_records(repo, "OrthoBench"),
              *[record(path) for path in sorted(Path(case["checkpoint"]).iterdir()) if path.is_file()]]
    phylo = Path(case["baseline_output"]) / "orthohmm_phylogeny"
    for directory in ("alignments", "gene_trees", "checkpoints", "species_tree_inference"):
        inputs.extend(record(path) for path in sorted((phylo / directory).rglob("*")) if path.is_file())
    inputs.extend(record(shutil.which(case[key]) or case[key]) for key in ("aligner", "tree_builder"))
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
            "repo": str(repo), "output": str(output), "cpu": 16, "arms": list(ARMS),
            "case": case, "inputs": inputs, "sources": source_pins, "runtime": scientific_runtime,
            "profile_runtime": require_profile_runtime(source), "tool_path": os.environ["PATH"],
            "scope": "Complete-proteome downstream HMM/profile inference from immutable canonical initial-search hits, with content-compatible tree-cache reuse. Shared-host incremental timing, not a fresh-search or isolated full-pipeline speed comparison.",
            "defaults_promoted": False}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", plan)
    return plan


def arm_commands(source, plan, folder, mode):
    if mode not in ARMS:
        raise ValueError("Unexpected weighting arm")
    case = plan["case"]
    graph, native = commands(source, plan, case, folder, Path(case["checkpoint"]), case["checkpoint_manifest"]["sha256"])
    if mode != "none":
        graph.extend(["--profile-sequence-weighting", mode])
    return graph, native


def validate_activation(graph, mode):
    parameters, counts = graph["parameters"], graph["counts"]
    if (parameters.get("profile_sequence_weighting") != mode
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
    if mode == "none":
        if effective:
            raise ValueError("Unweighted control unexpectedly records weighted profiles")
    elif (len(effective) != counts["profiles_built"]
          or len({item["cluster_id"] for item in effective}) != len(effective)
          or any(not math.isfinite(item["effective_sequence_count"])
                 or not 0 < item["effective_sequence_count"] <= 200 + 1e-8 for item in effective)):
        raise ValueError("Invalid per-profile weighting activation/effective counts")
    return {"mode": mode, "profiles_built": counts["profiles_built"],
            "weighted_profiles": len(effective), "strict_profile_edges": counts["strict_profile_edges"]}


def normalized(row):
    return {"f1": row["scores"]["OrthoBench F1"],
            "precision": row["details"]["precision"], "recall": row["details"]["recall"]}


def run(plan_path):
    plan = json.loads(plan_path.read_text())
    source, repo, output = (Path(plan[key]) for key in ("source_root", "repo", "output"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_profile_weighting_pilot.py":
        raise ValueError("Execute frozen pilot source")
    if tuple(plan["arms"]) != ARMS or plan["cpu"] != 16 or (output / "results.json").exists():
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
              "checkpoint_audit": checkpoint_audit, "ownership": ownership,
              "host_before": host_activity(), "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", report)
    control = None
    for mode in ARMS:
        row = {"dataset": "OrthoBench", "label": mode, "status": "running", "scores": {}}
        report["arms"].append(row)
        if mode != "none" and control is None:
            row["status"] = "not_run_failed_control_gate"
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
            row["activation"] = validate_activation(graph_report, mode)
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
            if mode == "none":
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
            if (native["genes"], native["species"]) != UNIVERSE:
                raise ValueError("Wrong native universe")
            row["native_summary"] = native["summary"]
            if mode == "none":
                row["golden_matches"] = compare_golden(phylo, case["golden"])
            row.update(score_orthobench(repo, prediction))
            actual, details, score_files = readback_scores(repo, "OrthoBench", row, native, folder, "weighting")
            if actual != row["scores"]:
                raise ValueError("Independent score readback differs")
            row.update(details=details, scoring_readback=score_files, metric_vector=endpoint_vector("OrthoBench", actual))
            for pin in (*checked, *native["checked"], *score_files):
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Runtime changed during evaluation")
            if mode == "none":
                if normalized(row) != case["expected"]:
                    raise ValueError("Replayed control scores differ from admitted canonical control")
                control = normalized(row)
            else:
                row["advancement"] = advancement(control, normalized(row))
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
    for name in ("repo", "source", "output"):
        preparation.add_argument("--" + name, required=True, type=Path)
    execution = commands_parser.add_parser("run")
    execution.add_argument("--execution-plan", required=True, type=Path)
    args = parser.parse_args(argv)
    result = (prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve())
              if args.command == "prepare" else run(args.execution_plan.resolve()))
    print(result["status"])
    return 0 if result["status"] in ("prepared_unrun", "finished_all_arms_verified") else 1


if __name__ == "__main__":
    raise SystemExit(main())

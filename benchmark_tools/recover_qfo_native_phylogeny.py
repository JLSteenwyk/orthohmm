"""Recover only native phylogeny/scoring after the original QfO child crash."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.report_hillclimb import readback_scores, replay_native
from benchmark_tools.run_full_search_hillclimb import checked_records, runtime
from benchmark_tools.run_held_tree_hillclimb import require_reuse
from benchmark_tools.run_membership_hillclimb import PRODUCTS, argument, score_qfo, validate_partition, write_json
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import sources, timed_stage, validate_ownership
from orthohmm.accuracy import load_accuracy_checkpoint

PROTOCOL = "hillclimb_qfo_native_recovery_20261014.json"
REVISION = "2ccc411520d9e63678178a50207ba1fca87bfbb3"


def require_native_crash(report, metrics):
    row = next(row for row in report["arms"] if row["label"] == "control100")
    if (report["status"] != "finished_with_explicit_gaps" or report["source_revision"] != REVISION
            or row["status"] != "failed" or row.get("inputs_sources_runtime_reverified") is not True
            or metrics.get("status") != "failed" or metrics.get("harness", {}).get("exit_code") != -11):
        raise ValueError("Require the authenticated terminal native signal-11 failure")


def prepare(repo, source, output):
    if output.exists() or source == repo:
        raise ValueError("Require new namespace and frozen recovery runner")
    parent_path = repo / "benchmarks/results/full_search_hillclimb_r1_v3_20261012/execution_plan.json"
    parent = json.loads(parent_path.read_text())
    if parent["source_revision"] != REVISION or runtime() != parent["runtime"]:
        raise ValueError("Original scientific source/runtime changed")
    case = next(case for case in parent["cases"] if case["dataset"] == "QfO" and case["label"] == "control100")
    original = Path(case["output"])
    report_path, metrics_path = original.parent / "results.json", original / "metrics.json"
    report, metrics = json.loads(report_path.read_text()), json.loads(metrics_path.read_text())
    require_native_crash(report, metrics)
    science = Path(parent["source_root"])
    checkpoint = original / "output/orthohmm_working_res/high_sensitivity_checkpoint"
    candidate = original / "output/orthohmm_working_res/phylogeny_candidate_superfamilies.txt"
    constraints = candidate.parent / "phylogeny_candidate_merges.json"
    load_membership_constraints(constraints, candidate)
    owners = {}
    for pin in case["inputs"]:
        for sequence in SeqIO.parse(pin["path"], "fasta"):
            if sequence.id in owners or not sequence.seq:
                raise ValueError("Duplicate or empty input gene")
            owners[sequence.id] = Path(pin["path"]).name
    names, species, *_ = load_accuracy_checkpoint(checkpoint, verify=True)
    ownership = validate_ownership(names, species, owners)
    if (ownership["genes"], ownership["species"]) != (984137, 78):
        raise ValueError("Wrong complete corrected QfO universe")
    read_partition(candidate, set(names))
    inputs = [record(source / "benchmark_tools" / PROTOCOL), record(parent_path), record(report_path), record(metrics_path),
              *checked_records(parent, case), record(candidate), record(constraints),
              *[record(p) for p in sorted(checkpoint.iterdir()) if p.is_file()]]
    phylo = original / "output/orthohmm_phylogeny"
    for directory in ("alignments", "gene_trees", "checkpoints", "species_tree_inference"):
        inputs.extend(record(p) for p in sorted((phylo / directory).rglob("*")) if p.is_file())
    source_pins = sources(source) + [record(source / "benchmark_tools/recover_qfo_native_phylogeny.slurm")]
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    for pin in source_pins:
        relative = Path(pin["path"]).relative_to(source)
        if relative.suffix == ".so": continue
        blob = subprocess.check_output(["git", "show", revision + ":" + str(relative)], cwd=repo)
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Recovery runner differs from committed source")
    command = [sys.executable, str(science / "benchmark_tools/replay_phylogeny.py"),
        "--fasta-directory", case["fasta"], "--candidate-clusters", str(candidate),
        "--membership-constraints", str(constraints), "--output-directory", str(output / "native"),
        "--json", str(output / "metrics.json"), "--cpu", "32", "--species-tree-mode", "infer",
        "--aligner", argument(case["command"], "--aligner"), "--tree-builder", argument(case["command"], "--tree-builder"),
        "--root-rule", "species_overlap", "--pair-rule", "positive_paralogy", "--species-tree-rooting", "min_variance",
        "--satellite-membership-rule", "high_confidence_pair", "--checkpoint-source", str(original / "output")]
    for pin in (*inputs, *source_pins): check(pin)
    result = {"status": "prepared_unrun_native_recovery", "repo": str(repo), "source_root": str(source),
              "source_revision": revision, "scientific_source": str(science), "scientific_revision": REVISION,
              "output": str(output), "command": command, "inputs": inputs, "sources": source_pins,
              "runtime": parent["runtime"], "tool_path": parent["tool_path"], "ownership": ownership,
              "checkpoint_audit": audit(checkpoint, record(checkpoint / "manifest.json")["sha256"]),
              "scope": "Cached native phylogeny recovery from fresh failed control; not a fresh uninterrupted full-pipeline timing."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def run(path):
    plan = json.loads(path.read_text())
    source, science, output, repo = (Path(plan[key]) for key in ("source_root", "scientific_source", "output", "repo"))
    if Path(__file__).resolve() != source / "benchmark_tools/recover_qfo_native_phylogeny.py" or (output / "results.json").exists():
        raise ValueError("Require frozen runner and preserve previous attempts")
    checked = [record(path), *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    result = {"status": "native_recovery_running", "execution_plan": record(path), "job_id": os.environ.get("SLURM_JOB_ID"),
              "scope": plan["scope"], "host_before": host_activity(), "scores": {}, "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", result)
    try:
        for pin in checked: check(pin)
        if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
            raise ValueError("Changed runtime or CUDA visibility")
        accounting = subprocess.check_output(["sacct", "-j", "24190_1", "-X", "-P", "-n", "-o", "JobID,State%30,ExitCode"], text=True)
        if "24190_1|FAILED|1:0" not in accounting:
            raise ValueError("Original task is not authoritatively terminal failed")
        result["parent_accounting"] = accounting
        env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(science), "PYTHONHASHSEED": "0",
               "CUDA_VISIBLE_DEVICES": "", "PYTHONFAULTHANDLER": "1", "OMP_NUM_THREADS": "1",
               "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "32"}
        result["execution"] = timed_stage(plan["command"], output, "native", science, env)
        if result["execution"]["returncode"]:
            raise RuntimeError("Standalone native phylogeny recovery failed; inspect fault-enabled log")
        phylo = output / "native/orthohmm_phylogeny"
        prediction = output / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
        candidate = Path(argument(plan["command"], "--candidate-clusters"))
        result["coverage"] = validate_partition(candidate, prediction)
        result["outputs"] = [record(phylo / name) for name in PRODUCTS]
        native = replay_native(result, {"command": plan["command"]})
        if (native["genes"], native["species"]) != (984137, 78):
            raise ValueError("Recovery changed the input universe")
        require_reuse(native["summary"])
        result["native_summary"] = native["summary"]
        result["native_checked"] = native["checked"]
        result["status"] = "native_recovered_scoring_pending"
        write_json(output / "results.json", result)
        scoring_id = "qnr_" + hashlib.sha256(str(output).encode()).hexdigest()[:12]
        result.update(score_qfo(repo, phylo, output, scoring_id))
        actual, details, files = readback_scores(repo, "QfO", result, native, output, "ohmm_hillclimb_r1_" + scoring_id)
        result.update(scores=actual, details=details, scoring_readback=files)
        for pin in (*checked, *native["checked"], *files): check(pin)
        if runtime() != plan["runtime"]:
            raise ValueError("Runtime changed during recovery")
        result["status"] = "cached_native_control_and_all_six_scores_verified"
    except Exception as error:
        result.update(status="failed_native_recovery", error_type=type(error).__name__, error=str(error))
    result["host_after"] = host_activity()
    write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    preparation = sub.add_parser("prepare")
    for name in ("repo", "source", "output"): preparation.add_argument("--" + name, type=Path, required=True)
    worker = sub.add_parser("run")
    worker.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve()) if args.command == "prepare" else run(args.execution_plan.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"] in {"prepared_unrun_native_recovery", "cached_native_control_and_all_six_scores_verified"} else 1)

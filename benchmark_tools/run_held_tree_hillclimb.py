"""Isolate the species-tree handoff with fixed broader QfO candidates and gene trees."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_species_tree_robustness import topology
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
from benchmark_tools.run_membership_hillclimb import (
    PRODUCTS, argument, compare_golden, score_qfo, validate_partition, write_json,
)
from benchmark_tools.run_search_recovery import host_activity, runtime_records
from benchmark_tools.run_three_kingdoms_membership_supplement import timed_stage

PLAN_NAME = "hillclimb_held_tree_20261012.json"
LABELS = ("broader_own_tree", "broader_baseline_tree")
DETERMINISTIC = frozenset(("VGNC F", "SwissTrees F", "TreeFam-A F", "EC", "GO"))


def tree_contrast(left, right):
    trees = [Phylo.read(path, "newick") for path in (left, right)]
    taxa = [sorted(t.name for t in tree.get_terminals()) for tree in trees]
    clades = [topology(tree) for tree in trees]
    if taxa[0] != taxa[1]:
        raise ValueError("Species-tree taxon universes differ")
    distance = len(clades[0] ^ clades[1])
    return {"taxa": len(taxa[0]), "rooted_clade_distance": distance,
            "normalized_rooted_clade_distance": distance / (len(clades[0]) + len(clades[1]))}


def supplied_command(original, source, folder, tree, cpu):
    argv = list(original)
    if (argument(argv, "--species-tree-mode") != "infer" or "--species-tree" in argv
            or argument(argv, "--satellite-membership-rule") != "high_confidence_pair"
            or "--unconstrained-membership" in argv):
        raise ValueError("Require unchanged inferred-tree direct-policy source")
    argument(argv, "--membership-constraints")
    checkpoint = argument(argv, "--output-directory")
    argv[:2] = [sys.executable, str(source / "benchmark_tools/replay_phylogeny.py")]
    for flag, value in (("--species-tree-mode", "supplied"), ("--cpu", cpu),
                        ("--output-directory", folder / "output"), ("--json", folder / "metrics.json"),
                        ("--checkpoint-source", checkpoint)):
        argument(argv, flag)
        argv[argv.index(flag) + 1] = str(value)
    argv.extend(["--species-tree", str(tree)])
    return argv


def require_reuse(summary):
    families, hits = summary["reconciled_families"], summary["checkpoint_hits"]
    if (isinstance(families, bool) or not isinstance(families, int) or families <= 0
            or isinstance(hits, bool) or not isinstance(hits, int) or hits != families):
        raise ValueError("Expected every reconciled raw gene tree to hit its checkpoint")


def require_control_scores(scores, expected):
    endpoint_vector("QfO", scores)
    endpoint_vector("QfO", expected)
    if any(scores[name] != expected[name] for name in DETERMINISTIC):
        raise ValueError("Deterministic own-tree control scores differ")


def source_records(source):
    return [record(path) for path in sorted((source / "orthohmm").rglob("*.py"))] + [
        record(path) for path in sorted((source / "benchmark_tools").glob("*.py"))] + [
        record(source / "qfo_benchmark/og_to_pairwise.py"),
        record(source / "benchmark_tools/run_held_tree_hillclimb.slurm")]


def prepare(repo, source, output, comparison, cpu):
    if output.exists() or isinstance(cpu, bool) or not 1 <= cpu <= 16:
        raise ValueError("Require fresh output and 1-16 CPUs")
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools"], check=True)
    report = json.loads(comparison.read_text())
    verified = {}
    for label in ("control_direct", "broader_direct"):
        method = next(m for m in report["methods"] if m["key"] == "membership:" + label)
        cells = [c for c in method["cells"] if c["dataset"] == "QfO"
                 and c["status"] == "independent_content_readback_verified"]
        if len(cells) != 1:
            raise ValueError("Need one independently verified QfO cell per source")
        verified[label] = cells[0]
    broad = verified["broader_direct"]
    producer = json.loads(Path(broad["source_snapshot"]["path"]).read_text())
    execution_record = producer["execution_plan"]
    check(execution_record)
    execution = json.loads(Path(execution_record["path"]).read_text())
    case = next(c for c in execution["cases"] if c["dataset"] == "QfO" and c["label"] == "broader_direct")
    checkpoint = Path(argument(case["command"], "--output-directory")) / "orthohmm_phylogeny"
    trees = {label: Path(cell["products"][PRODUCTS[0]]["path"]).parent / "species_tree.rooted.nwk"
             for label, cell in verified.items()}
    native = json.loads((checkpoint / "provenance_manifest.json").read_text())
    if any(native[key] != value for key, value in (("root_duplication_rule", "species_overlap"),
            ("pair_orthology_rule", "positive_paralogy"), ("species_tree_rooting", "min_variance"))):
        raise ValueError("Unexpected source reconciliation rules")
    contrast = tree_contrast(trees["control_direct"], trees["broader_direct"])
    if contrast["taxa"] != 78:
        raise ValueError("Expected all 78 corrected QfO species")
    inputs = [record(comparison), record(source / "benchmark_tools" / PLAN_NAME), execution_record,
              *[record(path) for path in trees.values()]]
    for cell in verified.values():
        inputs.extend([cell["source_snapshot"], *cell["checked"]])
    # Raw trees and checkpoints are input to the replay, not only its old summary.
    for name in ("alignments", "checkpoints", "gene_trees"):
        paths = sorted((checkpoint / name).rglob("*"))
        if not any(path.is_file() for path in paths):
            raise ValueError("Missing checkpoint input directory: " + name)
        inputs.extend(record(path) for path in paths if path.is_file())
    for flag in ("--aligner", "--tree-builder"):
        executable = shutil.which(argument(case["command"], flag))
        if executable is None:
            raise FileNotFoundError(argument(case["command"], flag))
        inputs.append(record(executable))
    inputs = list({item["path"]: item for item in inputs}.values())
    for item in inputs:
        check(item)
    nonce = hashlib.sha256(str(output).encode()).hexdigest()[:16]
    cases = [{"label": label, "command": supplied_command(case["command"], source, output / label, trees[tree_label], cpu),
              "output": str(output / label), "tree": record(trees[tree_label]),
              "scoring_id": f"heldtree_{nonce}_{index}"}
             for index, (label, tree_label) in enumerate(zip(LABELS, ("broader_direct", "control_direct")))]
    result = {"status": "prepared_unrun", "repo": str(repo), "source_root": str(source), "output": str(output),
        "source_revision": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "cpu": cpu, "inputs": inputs, "sources": source_records(source), "runtime": runtime_records(),
        "tool_path": os.environ["PATH"], "cases": cases, "tree_contrast": contrast,
        "golden": broad["products"], "expected_control_scores": {key.split(": ", 1)[1]: value for key, value in broad["scores"].items()},
        "baseline_scores": verified["control_direct"]["scores"],
        "scope": "Development-exposed fixed-candidate cached-gene-tree reconciliation; not full-pipeline timing or default promotion."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def run(plan_path):
    identity = record(plan_path)
    plan = json.loads(plan_path.read_text())
    source, repo, output = (Path(plan[key]) for key in ("source_root", "repo", "output"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_held_tree_hillclimb.py":
        raise ValueError("Execute the frozen source")
    if tuple(c["label"] for c in plan["cases"]) != LABELS or (output / "results.json").exists():
        raise ValueError("Unexpected case inventory or existing results")
    checked = [identity, *plan["sources"], *plan["inputs"], *plan["runtime"]["files"]]
    for item in checked:
        check(item)
    report = {"status": "running", "execution_plan": identity, "source_revision": plan["source_revision"],
        "arms": [], "job_id": os.environ.get("SLURM_JOB_ID"), "defaults_promoted": False,
        "scope": plan["scope"], "host_activity_before": host_activity()}
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "PYTHONHASHSEED": "0",
           "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    control_valid = False
    for case in plan["cases"]:
        row = {"label": case["label"], "status": "running", "scores": {}}
        report["arms"].append(row)
        write_json(output / "results.json", report)
        started = time.monotonic()
        try:
            if case["label"] != LABELS[0] and not control_valid:
                row.update(status="not_run_control_gate_failed")
                continue
            for item in checked:
                check(item)
            folder = Path(case["output"])
            folder.mkdir(exist_ok=False)
            write_json(folder / "preflight.json", case)
            row["native_execution"] = timed_stage(case["command"], folder, "inference", source, env)
            if row["native_execution"]["returncode"]:
                raise RuntimeError("Native replay failed")
            phylo = folder / "output/orthohmm_phylogeny"
            prediction = folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"
            row["coverage"] = validate_partition(Path(argument(case["command"], "--candidate-clusters")), prediction)
            row["outputs"] = [record(phylo / name) for name in PRODUCTS]
            native = replay_native(row, case)
            if native["genes"] != 984137 or native["species"] != 78:
                raise ValueError("Unexpected corrected QfO gene/species universe")
            require_reuse(native["summary"])
            row["tree_handoff"] = tree_contrast(case["tree"]["path"], phylo / "species_tree.rooted.nwk")
            if row["tree_handoff"]["rooted_clade_distance"] != 0:
                raise ValueError("Supplied rooted topology was not preserved")
            if case["label"] == LABELS[0]:
                row["golden_matches"] = compare_golden(phylo, plan["golden"])
            row["native_summary"] = native["summary"]
            row["native_readback"] = {key: native[key] for key in ("genes", "species", "groups", "pairs", "checked")}
            for item in [*checked, *native["checked"]]:
                check(item)
            row["status"] = "native_outputs_verified_scoring_pending"
            write_json(output / "results.json", report)
            row.update(score_qfo(repo, phylo, folder, case["scoring_id"]))
            scores, details, files = readback_scores(repo, "QfO", row, native, folder,
                                                    "ohmm_hillclimb_r1_" + case["scoring_id"])
            row.update(scores=scores, details=details, scoring_readback=files)
            if case["label"] == LABELS[0]:
                require_control_scores(scores, plan["expected_control_scores"])
            for item in [*checked, *native["checked"], *files]:
                check(item)
            row.update(status="independent_content_readback_verified", inputs_reverified_unchanged=True)
            if case["label"] == LABELS[0]:
                control_valid = True
        except Exception as error:
            row.update(status="failed", error_type=type(error).__name__, error=str(error))
        finally:
            row["total_incremental_wall_s"] = time.monotonic() - started
            try:
                for item in checked:
                    check(item)
            except Exception as error:
                row.update(status="failed", error_type=type(error).__name__, error=str(error))
                control_valid = False
            write_json(output / "results.json", report)
    report["status"] = ("finished_all_arms_verified" if all(row["status"] == "independent_content_readback_verified"
                                                           for row in report["arms"]) else "finished_with_explicit_gaps")
    report["host_activity_after"] = host_activity()
    write_json(output / "results.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="mode", required=True)
    preparation = commands.add_parser("prepare")
    for name in ("repo", "source-root", "output", "comparison"):
        preparation.add_argument("--" + name, type=Path, required=True)
    preparation.add_argument("--cpu", type=int, default=16)
    execution = commands.add_parser("run")
    execution.add_argument("--plan", type=Path, required=True)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source_root.resolve(), args.output.resolve(), args.comparison.resolve(), args.cpu)
              if args.mode == "prepare" else run(args.plan.resolve()))
    print(json.dumps({"status": result["status"]}))
    return int(args.mode == "run" and result["status"] != "finished_all_arms_verified")


if __name__ == "__main__":
    raise SystemExit(main())

"""Execute the prespecified fresh-search cap trial through native phylogeny."""

from __future__ import annotations

import argparse
import csv
import importlib
from itertools import zip_longest
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from orthohmm.accuracy import load_accuracy_checkpoint
from benchmark_tools.admit_qfo_factorial_cell import check_pairs
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import (
    DATASETS, PRODUCTS, argument, original_commands, scoring_records,
    score_orthobench, score_qfo, score_three_kingdoms, write_json,
)
from benchmark_tools.run_search_recovery import host_activity, runtime_records
from benchmark_tools.run_three_kingdoms_membership_supplement import validate_ownership
from benchmark_tools.score_ygob_groups import read_predictions, membership
from benchmark_tools.validate_profile_runtime import require_profile_runtime

PLAN_NAME = "hillclimb_full_search_20261012.json"


def runtime():
    result = runtime_records()
    for name in ("scipy", "igraph", "leidenalg"):
        module = importlib.import_module(name)
        result["versions"][name] = getattr(module, "__version__", "unspecified")
        result["files"].append(record(module.__file__))
    return result


def native_command(source, fasta, folder, tools, cap, cpu):
    return [sys.executable, str(source / "benchmark_tools/benchmark_production.py"),
        str(fasta), str(folder / "output"), str(folder / "metrics.json"),
        "--cpu", str(cpu), "--threads-per-worker", "8", "--full-output",
        "--matrix", "BLOSUM62", "--evalue", "0.0001", "--clustering", "leiden",
        "--cpm-resolution", "0.1", "--refinement-profile", "default",
        "--accuracy-profile", "high_sensitivity", "--search-max-candidates-per-query", str(cap),
        "--phylogeny", "reconcile", "--species-tree-mode", "infer",
        "--aligner", tools[0], "--tree-builder", tools[1],
        "--phylogeny-candidates", "satellite_v2", "--phylogeny-root-rule", "species_overlap",
        "--phylogeny-pair-rule", "positive_paralogy", "--species-tree-rooting", "min_variance"]


def prepare(repo, source, output, cpu=32):
    if output.exists() or cpu != 32:
        raise ValueError("Require fresh namespace and prespecified 32 CPUs")
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools"], check=True)
    plan_path = source / "benchmark_tools" / PLAN_NAME
    plan = json.loads(plan_path.read_text())
    if plan["arms"] != [{"label": "control100", "search_max_candidates_per_query": 100},
                        {"label": "cap400", "search_max_candidates_per_query": 400}]:
        raise ValueError("Unexpected prespecified arms")
    paths = [*sorted((source / "orthohmm").rglob("*.py")),
             *sorted((source / "orthohmm/search/csrc").glob("*.c")),
             *sorted((source / "orthohmm/search/csrc").glob("*.so")),
             *sorted((source / "benchmark_tools").glob("*.py")),
             source / "benchmark_tools/run_full_search_hillclimb.slurm"]
    for name in ("pair_align.so", "hmm_viterbi.so", "kmer_prefilter.so"):
        if not (source / "orthohmm/search/csrc" / name).is_file():
            raise ValueError("Missing CPU extension: " + name)
    cases = []
    for dataset in DATASETS:
        historical, _, evidence = original_commands(repo, dataset)
        fasta = Path(argument(historical, "--fasta-directory"))
        tools = [argument(historical, flag) for flag in ("--aligner", "--tree-builder")]
        tool_records = [record(shutil.which(tool) or tool) for tool in tools]
        inputs = [record(p) for p in sorted(fasta.iterdir()) if p.is_file()]
        expected_species = {"OrthoBench": 12, "QfO": 78, "ThreeKingdoms": 12}[dataset]
        if len(inputs) != expected_species:
            raise ValueError("Unexpected input proteome inventory: " + dataset)
        baseline = Path(argument(historical, "--output-directory")) / "orthohmm_phylogeny"
        golden = {name: record(baseline / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")}
        for index, arm in enumerate(plan["arms"]):
            folder = output / dataset / arm["label"]
            cases.append({"dataset": dataset, **arm, "fasta": str(fasta), "inputs": inputs,
                "evidence": evidence, "tools": tool_records, "golden": golden,
                "scorers": scoring_records(repo, dataset), "output": str(folder),
                "scoring_id": f"fullsearch_20261012_{dataset}_{index}",
                "command": native_command(source, fasta, folder, tools, arm["search_max_candidates_per_query"], cpu)})
    report = {"schema_version": 1, "status": "prepared_unrun", "repo": str(repo),
        "source_root": str(source), "source_revision": subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "output": str(output), "cpu": cpu, "plan": record(plan_path),
        "selection_evidence": record(repo / plan["selection_evidence"]),
        "sources": [record(p) for p in paths], "runtime": runtime(), "cases": cases,
        "tool_path": os.environ["PATH"], "scope": "Fresh full-proteome CPU HMM search through inferred phylogeny; shared-host development comparison."}
    for item in checked_records(report):
        check(item)
    report["profile_runtime"] = require_profile_runtime(source)
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", report)
    return report


def checked_records(plan, case=None):
    records = [plan["plan"], plan["selection_evidence"], *plan["sources"], *plan["runtime"]["files"]]
    for row in ([case] if case is not None else plan["cases"]):
        records.extend([*row["inputs"], *row["evidence"], *row["tools"], *row["golden"].values(), *row["scorers"]])
    return records


def validate_confidence(path, pairs):
    with pairs.open() as left, path.open() as right:
        a, b = csv.reader(left, delimiter="\t"), csv.reader(right, delimiter="\t")
        header = ["gene_a", "species_a", "gene_b", "species_b"]
        if next(a, None) != header or next(b, None) != [*header, "confidence"]:
            raise ValueError("Invalid pair/confidence header")
        count = 0
        for pair, confidence in zip_longest(a, b):
            if (pair is None or confidence is None or len(confidence) != 5
                    or confidence[:4] != pair or confidence[4] not in {"high", "medium", "low"}):
                raise ValueError("Pair/confidence outputs differ")
            count += 1
    return count


def validate_native(case, cpu):
    folder = Path(case["output"])
    output, phylo = folder / "output", folder / "output/orthohmm_phylogeny"
    metrics = json.loads((folder / "metrics.json").read_text())
    if metrics.get("status") != "complete" or metrics.get("harness", {}).get("exit_code") != 0:
        raise ValueError("Incomplete fresh production run")
    metadata = metrics["metadata"]
    expected = {"accuracy_profile": "high_sensitivity", "search_mode": "builtin", "clustering": "leiden",
        "cpm_resolution": .1, "substitution_matrix": "BLOSUM62", "evalue_threshold": 1e-4,
        "leiden_seed": 4, "cpu_budget": cpu, "search_kmer_k": 4,
        "search_max_candidates_per_query": case["search_max_candidates_per_query"],
        "phylogeny": "reconcile", "phylogeny_candidates": "satellite_v2", "species_tree_mode": "infer",
        "phylogeny_root_rule": "species_overlap", "phylogeny_pair_rule": "positive_paralogy",
        "species_tree_rooting": "min_variance", "satellite_membership_rule": "high_confidence_pair"}
    if any(metadata.get(key) != value for key, value in expected.items()):
        raise ValueError("Actual native parameters differ from trial")
    if Path(metadata["fasta_directory"]).resolve() != Path(case["fasta"]).resolve() or Path(metadata["output_directory"]).resolve() != output.resolve():
        raise ValueError("Wrong input/output directory metadata")
    native = json.loads((phylo / "provenance_manifest.json").read_text())
    summary = json.loads((phylo / "reconciliation_summary.json").read_text())
    taxa = {item["filename"]: item["taxon"] for item in native["input_proteomes"]}
    if len(taxa) != len(case["inputs"]) or len(set(taxa.values())) != len(taxa):
        raise ValueError("Nonbijective native species inventory")
    owners, filenames = {}, {}
    for item in case["inputs"]:
        path = Path(item["path"])
        if path.name not in taxa:
            raise ValueError("Unexpected native/input FASTA inventory")
        for seq in SeqIO.parse(path, "fasta"):
            if seq.id in owners or not seq.seq:
                raise ValueError("Duplicate or empty input sequence")
            owners[seq.id], filenames[seq.id] = taxa[path.name], path.name
    checkpoint = output / "orthohmm_working_res/high_sensitivity_checkpoint"
    numeric = audit(checkpoint, record(checkpoint / "manifest.json")["sha256"])
    names, species, *_ = load_accuracy_checkpoint(checkpoint, verify=False)
    ownership = validate_ownership(names, species, filenames)
    candidate_path = output / "orthohmm_working_res/phylogeny_candidate_superfamilies.txt"
    candidates = {g: i for i, group in enumerate(read_partition(candidate_path, set(owners))) for g in group}
    groups = read_partition(output / "orthohmm_working_res/orthohmm_edges_clustered.txt", set(owners))
    root = read_predictions(phylo / "orthohmm_root_hogs.tsv", "root_hogs")
    exported = read_predictions(output / "orthohmm_orthogroups.txt", "named_groups")
    if (set(membership(root)) != set(owners) or set(membership(exported)) != set(owners)
            or {frozenset(g) for g in groups} != {frozenset(g) for g in root.values()}
            or {frozenset(g) for g in groups} != {frozenset(g) for g in exported.values()}):
        raise ValueError("Native partition/root/exported products differ")
    count = check_pairs(phylo / PRODUCTS[1], owners, candidates, summary["ortholog_pairs"])
    if count <= 0 or validate_confidence(phylo / PRODUCTS[2], phylo / PRODUCTS[1]) != count:
        raise ValueError("Invalid native pair/confidence count")
    if (summary["checkpoint_hits"] != 0 or summary["species_tree_checkpoint_hit"] is not False
            or summary.get("remapped_checkpoint_hits", 0) != 0 or summary["reconciled_families"] <= 0
            or summary["species_tree_families"] <= 0):
        raise ValueError("Require fresh inferred trees and reconciliation")
    counts = metrics["counts"]
    expected_counts = {"genes": len(owners), "species": len(taxa), "species_pairs": len(taxa)**2,
                       "orthogroups": len(groups), "phylogeny_ortholog_pairs": count}
    if any(counts.get(key) != value for key, value in expected_counts.items()):
        raise ValueError("Native counts differ from independently read outputs")
    if any(type(counts.get(key)) is not int or counts[key] <= 0 for key in
           ("search_candidates", "significant_hits", "high_sensitivity_profiles")):
        raise ValueError("Missing initial/profile HMM evidence")
    products = {name: record(phylo / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")}
    return {"status": "fresh_native_content_verified", "genes": len(owners), "species": len(taxa),
        "groups": len(groups), "pairs": count, "numeric_checkpoint": numeric, "ownership": ownership,
        "historical_product_matches": {name: item["sha256"] == case["golden"][name]["sha256"] for name, item in products.items()},
        "products": products, "metrics": metrics, "summary": summary}


def run(path, dataset):
    plan = json.loads(path.read_text())
    source, repo = Path(plan["source_root"]), Path(plan["repo"])
    if Path(__file__).resolve() != source / "benchmark_tools/run_full_search_hillclimb.py":
        raise ValueError("Execute frozen runner")
    folder = Path(plan["output"]) / dataset
    folder.mkdir(exist_ok=False)
    report = {"status": "running", "dataset": dataset, "execution_plan": record(path),
        "source_revision": plan["source_revision"], "job_id": os.environ.get("SLURM_JOB_ID"),
        "arms": [], "defaults_promoted": False, "scope": plan["scope"]}
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "PYTHONHASHSEED": "0",
        "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        "NUMBA_NUM_THREADS": str(plan["cpu"])}
    control_valid = False
    for case in [c for c in plan["cases"] if c["dataset"] == dataset]:
        row = {"label": case["label"], "status": "running", "scores": {}}
        report["arms"].append(row)
        write_json(folder / "results.json", report)
        started = time.monotonic()
        try:
            for item in [record(path), *checked_records(plan, case)]:
                check(item)
            if runtime() != plan["runtime"]:
                raise ValueError("Runtime changed")
            if case["label"] != "control100" and not control_valid:
                row["status"] = "not_run_failed_control"
                continue
            arm = Path(case["output"])
            arm.mkdir(exist_ok=False)
            (arm / "host_before.txt").write_text(host_activity())
            row["profile_runtime"] = require_profile_runtime(source)
            probe = subprocess.check_output([sys.executable, "-c",
                "from orthohmm.search import engine; from orthohmm.search.prefilter import is_c_prefilter_available; "
                "assert not engine.is_cuda_available() and engine.is_c_available() and is_c_prefilter_available()"], cwd=source, env=env, text=True)
            row["cpu_backend_probe"] = {"passed": True, "stdout": probe}
            with (arm / "inference.log").open("x") as log:
                child = subprocess.run(case["command"], cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT)
            row["exit_code"] = child.returncode
            if child.returncode:
                raise RuntimeError("Fresh pipeline failed: " + str(child.returncode))
            row["native"] = validate_native(case, plan["cpu"])
            row["status"] = "native_verified_scoring_pending"
            write_json(folder / "results.json", report)
            prediction = arm / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"
            scoring = (score_orthobench(repo, prediction) if dataset == "OrthoBench" else
                       score_three_kingdoms(repo, prediction, arm) if dataset == "ThreeKingdoms" else
                       score_qfo(repo, arm / "output/orthohmm_phylogeny", arm, case["scoring_id"]))
            row.update(status="scored_development_result", **scoring)
            if case["label"] == "control100":
                control_valid = True
        except Exception as error:
            row.update(status="failed", error_type=type(error).__name__, error=str(error))
            control_valid = False
        finally:
            row["total_full_pipeline_and_scoring_wall_s"] = time.monotonic() - started
            try:
                for item in [report["execution_plan"], *checked_records(plan, case)]:
                    check(item)
                if runtime() != plan["runtime"]:
                    raise ValueError("Runtime changed during run")
                row["inputs_sources_runtime_reverified"] = True
            except Exception as error:
                row.update(status="failed", error_type=type(error).__name__, error=str(error))
                control_valid = False
            arm = Path(case["output"])
            if arm.is_dir():
                (arm / "host_after.txt").write_text(host_activity())
            write_json(folder / "results.json", report)
    report["status"] = "finished_all_arms_scored" if len(report["arms"]) == 2 and all(
        r["status"] == "scored_development_result" for r in report["arms"]) else "finished_with_explicit_gaps"
    write_json(folder / "results.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    prep = modes.add_parser("prepare")
    for name in ("repo", "source-root", "output"):
        prep.add_argument("--" + name, type=Path, required=True)
    worker = modes.add_parser("run")
    worker.add_argument("--execution-plan", type=Path, required=True)
    worker.add_argument("--dataset", choices=DATASETS, required=True)
    args = parser.parse_args()
    if args.mode == "prepare":
        result = prepare(args.repo.resolve(), args.source_root.resolve(), args.output.resolve())
    else:
        result = run(args.execution_plan.resolve(), args.dataset)
    print(json.dumps({"status": result["status"]}))
    return 0 if result["status"] in {"prepared_unrun", "finished_all_arms_scored"} else 1


if __name__ == "__main__":
    raise SystemExit(main())

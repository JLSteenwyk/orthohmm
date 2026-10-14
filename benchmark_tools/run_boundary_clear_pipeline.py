"""Run private boundary-clear libraries through matched fresh native inference."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.report_hillclimb import endpoint_vector
from benchmark_tools.run_full_search_hillclimb import (
    native_command, runtime, validate_native, verify_input_inventory,
)
from benchmark_tools.run_membership_hillclimb import argument, score_orthobench, scoring_records, write_json
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import timed_stage
from benchmark_tools.validate_profile_runtime import require_profile_runtime

ARMS = ("full_clear", "boundary_only")
PLAN = "hillclimb_boundary_clear_pipeline_20261014.json"
ARTIFACTS = (
    "orthohmm_working_res/high_sensitivity_checkpoint/manifest.json",
    "orthohmm_working_res/orthohmm_edges.txt",
    "orthohmm_working_res/orthohmm_edges_clustered.txt",
    "orthohmm_working_res/phylogeny_candidate_superfamilies.txt",
    "orthohmm_working_res/phylogeny_candidate_merges.json",
    "orthohmm_working_res/phylogeny_candidate_seeds.tsv",
    "orthohmm_phylogeny/orthohmm_root_hogs.tsv",
    "orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv",
    "orthohmm_phylogeny/orthohmm_pairwise_orthologs_confidence.tsv",
    "orthohmm_phylogeny/species_tree.rooted.nwk",
    "orthohmm_orthogroups.txt",
    "orthohmm_gene_count.txt",
)


def route_totals(directory, expected_library, counts):
    files = sorted(directory.glob("routes_*.json"))
    totals = {route: {"calls": 0, "scored_pairs": 0} for route in ("initial_search", "profile_scoring")}
    pids = set()
    for path in files:
        report = json.loads(path.read_text())
        if report["pid"] in pids or path.name != f"routes_{report['pid']}.json" or report["library"] != expected_library:
            raise ValueError("Invalid or mixed actual native route evidence")
        pids.add(report["pid"])
        for name, stats in report["routes"].items():
            if (name not in totals or set(stats) != {"calls", "scored_pairs"}
                    or type(stats["calls"]) is not int or stats["calls"] <= 0
                    or type(stats["scored_pairs"]) is not int or stats["scored_pairs"] < 0):
                raise ValueError("Invalid route counts")
            for key in stats:
                totals[name][key] += stats[key]
    if (totals["initial_search"]["scored_pairs"] != counts["search_candidates"]
            or totals["profile_scoring"]["scored_pairs"] < counts["high_sensitivity_profile_candidates"]
            or any(row["calls"] <= 0 for row in totals.values())):
        raise ValueError("Actual native route counts do not cover initial/profile inference")
    return {"totals": totals, "files": [record(path) for path in files]}


def product_records(output):
    records = {name: record(output / name) for name in ARTIFACTS}
    for name in ("orthohmm_orthogroups", "orthohmm_single_copy_orthogroups"):
        directory = output / name
        if not directory.is_dir():
            raise ValueError("Missing full-output sequence export directory")
        records.update({str(path.relative_to(output)): record(path) for path in sorted(directory.rglob("*")) if path.is_file()})
    return records


def compare_products(left, right):
    if set(left) != set(right):
        raise ValueError("Full biological export inventories differ")
    changed = [name for name in left if (left[name]["sha256"], left[name]["bytes"]) != (right[name]["sha256"], right[name]["bytes"])]
    if changed:
        raise ValueError("Full native output equivalence failed: " + str(changed[:20]))
    return {"files": len(left), "all_hashes_equal": True}


def prepare(repo, source, output):
    if output.exists() or source == repo or (source / "variants").exists():
        raise ValueError("Require fresh output and private frozen source namespace")
    decision_path = repo / "benchmark_tools/results/boundary_clear_search_terminal_decision_20261014.json"
    decision = json.loads(decision_path.read_text())
    if decision["status"] != "broader_complete_search_output_equivalence_and_lower_observed_cost_verified":
        raise ValueError("Complete-search kernel gate not admitted")
    c_path = source / "orthohmm/search/csrc/hmm_viterbi.c"
    previous = subprocess.check_output(["git", "show", decision["source_revision"] + ":orthohmm/search/csrc/hmm_viterbi.c"], cwd=repo)
    if hashlib.sha256(previous).hexdigest() != record(c_path)["sha256"]:
        raise ValueError("Kernel changed after complete-search admission")
    parent_path = repo / "benchmarks/results/full_search_hillclimb_r1_v3_20261012/execution_plan.json"
    parent = json.loads(parent_path.read_text())
    old = next(case for case in parent["cases"] if case["dataset"] == "OrthoBench" and case["label"] == "control100")
    golden = old["golden"]
    inputs = [record(source / "benchmark_tools" / PLAN), record(decision_path), record(parent_path),
              *old["inputs"], *old["tools"], *golden.values(), *scoring_records(repo, "OrthoBench")]
    compiler = shutil.which("gcc")
    if compiler is None:
        raise FileNotFoundError("GCC unavailable")
    inputs.append(record(compiler))
    original = Path(parent["source_root"])
    for name in ("kmer_prefilter.so", "pair_align.so"):
        inputs.append(record(original / "orthohmm/search/csrc" / name))
    tool_names = [argument(old["command"], flag) for flag in ("--aligner", "--tree-builder")]
    cases, source_pins = [], []
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    for arm in ARMS:
        variant = source / "variants" / arm
        shutil.copytree(source / "orthohmm", variant / "orthohmm", ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.so"))
        (variant / "benchmark_tools").mkdir()
        for path in (source / "benchmark_tools").glob("*.py"):
            shutil.copy2(path, variant / "benchmark_tools" / path.name)
        for name in ("kmer_prefilter.so", "pair_align.so"):
            shutil.copy2(original / "orthohmm/search/csrc" / name, variant / "orthohmm/search/csrc" / name)
        kernel = variant / "orthohmm/search/csrc/hmm_viterbi.so"
        flag = ["-DHMM_VITERBI_BOUNDARY_CLEAR_ONLY=1"] if arm == "boundary_only" else []
        command = [compiler, "-O3", "-fopenmp", "-shared", "-fPIC", "-march=native", *flag,
                   str(variant / "orthohmm/search/csrc/hmm_viterbi.c"), "-o", str(kernel)]
        subprocess.run(command, check=True, capture_output=True)
        folder = output / arm
        native = native_command(variant, Path(old["fasta"]), folder, tool_names, 100, 16)
        native.extend(["--search-result-order", "canonical", "--native-kernel-trace", str(folder / "native_routes")])
        cases.append({"label": arm, "source_root": str(variant), "output": str(folder), "fasta": old["fasta"],
                      "inputs": old["inputs"], "golden": golden, "search_max_candidates_per_query": 100,
                      "kernel": record(kernel), "compiler_command": command, "command": native,
                      "profile_runtime": require_profile_runtime(variant)})
        source_pins.extend(record(path) for path in sorted((variant / "orthohmm").rglob("*")) if path.is_file() and path.suffix in {".py", ".c", ".so"})
        source_pins.extend(record(path) for path in sorted((variant / "benchmark_tools").glob("*.py")))
    base_sources = [record(path) for path in sorted((source / "benchmark_tools").glob("*.py"))]
    base_sources.extend(record(path) for path in sorted((source / "orthohmm").rglob("*.py")))
    base_sources.append(record(source / "benchmark_tools/run_boundary_clear_pipeline.slurm"))
    for pin in base_sources:
        relative = Path(pin["path"]).relative_to(source)
        blob = subprocess.check_output(["git", "show", revision + ":" + str(relative)], cwd=repo)
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Frozen runner source changed")
    # Variant scientific sources are copies of the admitted frozen source, not edits.
    for pin in source_pins:
        path = Path(pin["path"])
        arm = next(case for case in cases if path.is_relative_to(Path(case["source_root"])))
        if path.suffix != ".so" and pin["sha256"] != record(source / path.relative_to(arm["source_root"]))["sha256"]:
            raise ValueError("Variant source differs beyond the compile flag")
    scientific_runtime = runtime()
    for pin in (*inputs, *source_pins, *base_sources, *scientific_runtime["files"]):
        check(pin)
    result = {"status": "prepared_unrun", "repo": str(repo), "source_root": str(source), "source_revision": revision,
              "output": str(output), "cpu": 16, "inputs": inputs, "sources": [*base_sources, *source_pins],
              "runtime": scientific_runtime, "tool_path": os.environ["PATH"], "cases": cases,
              "compiler_version": subprocess.check_output([compiler, "--version"], text=True),
              "scope": "Fresh full native OrthoBench output equivalence; one ordered shared-host pair is descriptive cost only.",
              "defaults_promoted": False}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def run(plan_path):
    plan = json.loads(plan_path.read_text())
    output, source, repo = (Path(plan[key]) for key in ("output", "source_root", "repo"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_boundary_clear_pipeline.py" or (output / "results.json").exists():
        raise ValueError("Require frozen runner and preserve previous attempts")
    checked = [record(plan_path), *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    result = {"status": "running", "job_id": os.environ.get("SLURM_JOB_ID"), "execution_plan": record(plan_path),
              "source_revision": plan["source_revision"], "scope": plan["scope"], "arms": [],
              "host_before": host_activity(), "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", result)
    control = None
    for case in plan["cases"]:
        row = {"label": case["label"], "status": "running", "scores": {}}
        result["arms"].append(row)
        write_json(output / "results.json", result)
        if case["label"] != ARMS[0] and control is None:
            row["status"] = "not_run_failed_control"
            write_json(output / "results.json", result)
            continue
        folder, variant = Path(case["output"]), Path(case["source_root"])
        folder.mkdir(exist_ok=False)
        try:
            for pin in checked:
                check(pin)
            if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
                raise ValueError("Changed runtime or CUDA visibility")
            verify_input_inventory(case)
            env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(variant), "PYTHONHASHSEED": "0",
                   "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "16"}
            row["execution"] = timed_stage(case["command"], folder, "native", variant, env)
            if row["execution"]["returncode"]:
                raise RuntimeError("Fresh strict native pipeline failed")
            row["native"] = validate_native(case, 16)
            if (row["native"]["genes"], row["native"]["species"]) != (251378, 12):
                raise ValueError("Wrong full OrthoBench universe")
            if row["native"]["metrics"]["metadata"].get("search_result_order") != "canonical":
                raise ValueError("Canonical collection not active")
            row["route_evidence"] = route_totals(folder / "native_routes", case["kernel"], row["native"]["metrics"]["counts"])
            row["products"] = product_records(folder / "output")
            row.update(score_orthobench(repo, folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"))
            row["metric_vector"] = endpoint_vector("OrthoBench", row["scores"])
            if control is not None:
                row["exact_product_comparison"] = compare_products(control["products"], row["products"])
                if control["native"]["metrics"]["counts"] != row["native"]["metrics"]["counts"]:
                    raise ValueError("Algorithmic pipeline counts differ")
                left = {k: v for k, v in control["native"]["summary"].items() if type(v) in (int, bool)}
                right = {k: v for k, v in row["native"]["summary"].items() if type(v) in (int, bool)}
                if left != right or row["scores"] != control["scores"] or row["details"] != control["details"]:
                    raise ValueError("Native reconciliation counts or independent scores differ")
            for pin in (*checked, *row["products"].values(), *row["route_evidence"]["files"], *row["references"]):
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Runtime changed during inference")
            row["status"] = "fresh_native_and_scores_verified"
            if case["label"] == ARMS[0]:
                control = row
        except Exception as error:
            row.update(status="failed", error_type=type(error).__name__, error=str(error))
            control = None
        write_json(output / "results.json", result)
        print(case["label"], row["status"], flush=True)
    result["status"] = ("complete_fresh_full_native_output_equivalence" if len(result["arms"]) == 2 and all(
        row["status"] == "fresh_native_and_scores_verified" for row in result["arms"]) else "finished_with_explicit_gaps")
    result["host_after"] = host_activity()
    write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    preparation = sub.add_parser("prepare")
    for name in ("repo", "source", "output"):
        preparation.add_argument("--" + name, required=True, type=Path)
    worker = sub.add_parser("run")
    worker.add_argument("--execution-plan", required=True, type=Path)
    args = parser.parse_args()
    result = prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve()) if args.command == "prepare" else run(args.execution_plan.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"] in {"prepared_unrun", "complete_fresh_full_native_output_equivalence"} else 1)

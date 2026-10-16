"""Finish the original cap400 trial from an authenticated recovered final graph."""

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_cap400_final_graph import authenticate, REVISION
from benchmark_tools.recover_full_qfo_search import short_id
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import PRODUCTS, argument, score_qfo, validate_partition, write_json
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import timed_stage

PROTOCOL = "hillclimb_cap400_downstream_recovery_20261016.json"
ADMISSION = "benchmark_tools/results/cap400_final_graph_partition_readback_20261016.json"


def require_graph_success(admission, plan, report, scheduler):
    if (admission.get("status") != "original_frozen_final_partition_independently_read_back_downstream_pending"
            or report.get("status") != "final_graph_partition_verified_downstream_pending"
            or report.get("optimizer_attempted") is not True
            or report.get("postflight_pins_verified") is not True
            or report.get("execution", {}).get("returncode") != 0
            or report.get("scores") != {} or admission.get("scores") != {}
            or admission.get("scientific_revision") != REVISION
            or plan.get("scientific_revision") != REVISION):
        raise ValueError("Require verified successful recovered graph, without accuracy claims")
    if (scheduler.get("JobIDRaw"), scheduler.get("State"), scheduler.get("ExitCode"),
            scheduler.get("AllocCPUS"), scheduler.get("NodeList")) != (
            "24347", "COMPLETED", "0:0", "1", "bizon"):
        raise ValueError("Graph recovery scheduler is not terminal successful")
    if (report.get("independent_rows") != {
            "status": "all_saved_rows_and_float64_bits_verified", "edges": 33585180,
            "genes": 984137, "isolates_preserved_by_complete_gene_table": True,
            "edge_rows_deduplicated_or_reordered": False}
            or admission.get("partition") != report.get("partition", {}).get("product")
            or admission.get("genes") != 984137 or admission.get("groups") != 276617):
        raise ValueError("Recovered graph row verification or partition admission differs")
    if plan.get("metadata") != {"cpm_resolution": .1, "seed": 4, "include_isolates": True,
                                "output_directory": str(Path(plan["output"]) / "native")}:
        raise ValueError("Recovered optimizer scientific settings differ")


def verify_graph(repo, admission_pin):
    check(admission_pin)
    admission = json.loads(Path(admission_pin["path"]).read_text())
    for pin in (admission["execution_plan"], admission["producer"], admission["partition"],
                admission["complete_gene_order"], admission["saved_graph"], *admission["reconstructed_arrays"]):
        check(pin)
    graph_plan = json.loads(Path(admission["execution_plan"]["path"]).read_text())
    graph_report = json.loads(Path(admission["producer"]["path"]).read_text())
    rows = list(csv.DictReader(io.StringIO(subprocess.check_output([
        "sacct", "-j", "24347", "-X", "-P", "-o", "JobIDRaw,State,ExitCode,AllocCPUS,NodeList"],
        text=True)), delimiter="|"))
    if len(rows) != 1:
        raise ValueError("Graph recovery scheduler record is not unique")
    require_graph_success(admission, graph_plan, graph_report, rows[0])
    if graph_report["execution"]["command"] != [sys.executable, "-m", "orthohmm.leiden_worker",
                                                str(Path(graph_plan["output"]) / "payload")]:
        raise ValueError("Recovered graph did not use original frozen worker command")
    parent, failure, _, _ = authenticate(repo, graph_plan["failure"])
    if (graph_report["execution_plan"] != admission["execution_plan"]
            or graph_plan["abc"] != admission["saved_graph"]
            or graph_report["reconstruction"]["files"] != admission["reconstructed_arrays"]):
        raise ValueError("Graph recovery provenance differs from independent admission")
    names = Path(admission["complete_gene_order"]["path"]).read_text().splitlines()
    if len(names) != 984137 or len(read_partition(Path(admission["partition"]["path"]), set(names))) != 276617:
        raise ValueError("Recovered partition universe or counts differ")
    return parent, failure, admission, graph_plan


def phylogeny_command(plan):
    working = Path(plan["output"]) / "graph/orthohmm_working_res"
    return [sys.executable, str(Path(plan["scientific_source"]) / "benchmark_tools/replay_phylogeny.py"),
        "--fasta-directory", plan["fasta"], "--candidate-clusters", str(working / "phylogeny_candidate_superfamilies.txt"),
        "--membership-constraints", str(working / "phylogeny_candidate_merges.json"),
        "--output-directory", str(Path(plan["output"]) / "native"),
        "--json", str(Path(plan["output"]) / "native_metrics.json"), "--cpu", "32",
        "--species-tree-mode", "infer", "--aligner", plan["aligner"], "--tree-builder", plan["tree_builder"],
        "--root-rule", "species_overlap", "--pair-rule", "positive_paralogy",
        "--species-tree-rooting", "min_variance", "--satellite-membership-rule", "high_confidence_pair"]


def require_original_settings(plan, parent, admission, graph_plan, original_result):
    case = parent["case"]["fresh_case"]
    expected = {"scientific_source": parent["scientific_source"], "scientific_revision": REVISION,
        "fasta": case["fasta"], "aligner": argument(case["command"], "--aligner"),
        "tree_builder": argument(case["command"], "--tree-builder"),
        "checkpoint": str(Path(admission["complete_gene_order"]["path"]).parent),
        "partition": admission["partition"], "payload": str(Path(graph_plan["output"]) / "payload"),
        "runtime": parent["runtime"], "tool_path": parent["tool_path"], "control": original_result["control"]}
    if any(plan.get(key) != value for key, value in expected.items()):
        raise ValueError("Downstream settings differ from authenticated original trial")
    if plan["command"] != phylogeny_command(plan):
        raise ValueError("Downstream command differs from original phylogeny configuration")


def prepare(repo, source, output):
    if source == repo or output.exists():
        raise ValueError("Require a frozen runner and new downstream namespace")
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    owned = [source / "benchmark_tools" / name for name in (PROTOCOL,
        "continue_cap400_downstream.py", "continue_cap400_downstream_worker.py", "continue_cap400_downstream.slurm")]
    for path in owned:
        blob = subprocess.check_output(["git", "-C", str(source), "show", revision + ":" + str(path.relative_to(source))])
        if hashlib.sha256(blob).hexdigest() != record(path)["sha256"]:
            raise ValueError("Continuation differs from committed source")
    admission_pin = record(repo / ADMISSION)
    if record(source / ADMISSION)["sha256"] != admission_pin["sha256"]:
        raise ValueError("Partition admission differs from frozen evidence")
    parent, failure, admission, graph_plan = verify_graph(repo, admission_pin)
    original_result = json.loads(Path(failure["parent_result"]["path"]).read_text())
    control = original_result["control"]
    endpoint_vector("QfO", control["scores"])
    checked = [admission_pin, failure["execution_plan"], failure["parent_result"],
        admission["execution_plan"], admission["producer"], admission["partition"], admission["saved_graph"],
        *admission["reconstructed_arrays"], *parent["inputs"], *parent["sources"], *parent["runtime"]["files"],
        control["report"], *control["checked"]]
    checkpoint = Path(admission["complete_gene_order"]["path"]).parent
    checked += [record(p) for p in sorted(checkpoint.iterdir())]
    sources = [record(p) for p in sorted((source / "benchmark_tools").glob("*.py"))]
    sources += [record(p) for p in sorted((source / "orthohmm").rglob("*.py"))]
    sources += [record(p) for p in owned]
    for pin in checked + sources:
        check(pin)
    if runtime() != parent["runtime"]:
        raise ValueError("Original scientific runtime changed")
    case = parent["case"]["fresh_case"]
    plan = {"status": "prepared_unrun_cap400_downstream", "repo": str(repo), "source_root": str(source),
        "source_revision": revision, "output": str(output), "scientific_source": parent["scientific_source"],
        "scientific_revision": REVISION, "admission": admission_pin, "checked": checked, "sources": sources,
        "runtime": parent["runtime"], "tool_path": parent["tool_path"], "control": control,
        "checkpoint": str(checkpoint), "partition": admission["partition"],
        "payload": str(Path(graph_plan["output"]) / "payload"), "fasta": case["fasta"],
        "aligner": argument(case["command"], "--aligner"), "tree_builder": argument(case["command"], "--tree-builder"),
        "scoring_id": short_id(repo, output, 1),
        "scope": "Cached final-graph continuation of original cap400; no repeated HMM or clustering; not uninterrupted timing."}
    plan["command"] = phylogeny_command(plan)
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", plan)
    return plan


def validate_graph_stage(stage, plan, plan_pin):
    if (stage.get("status") != "original_refinement_and_satellite_candidates_verified"
            or stage.get("execution_plan") != plan_pin
            or stage.get("scientific_revision") != REVISION
            or (stage.get("ownership", {}).get("genes"), stage.get("ownership", {}).get("species")) != (984137, 78)
            or stage.get("initial_hmm_search_rerun") is not False
            or stage.get("profile_expansion_rerun") is not False
            or stage.get("graph_clustering_rerun") is not False or stage.get("scores") != {}):
        raise ValueError("Incomplete or changed original refinement/candidate execution")
    expected = Path(plan["output"]) / "graph/orthohmm_working_res/phylogeny_candidate_superfamilies.txt"
    if (stage["candidate"].get("profile") != "satellite_v2"
            or stage["candidate"].get("candidate_checkpoint") != str(expected)
            or stage["candidate"].get("membership_policy") != "high_confidence_pair"
            or stage["candidate"].get("seed_families") != stage.get("refined_groups")
            or stage["candidate"].get("candidate_families") != stage.get("candidate_groups")):
        raise ValueError("Changed native candidate profile or membership policy")
    graph = Path(plan["output"]) / "graph"
    expected_products = [graph / "orthogroups_refined.txt", expected,
        expected.parent / "phylogeny_candidate_merges.json", expected.parent / "phylogeny_candidate_seeds.tsv",
        expected.parent / "orthohmm_edges_clustered.txt"]
    expected_modules = [Path(plan["scientific_source"]) / "orthohmm" / name
                        for name in ("orthohmm.py", "accuracy.py", "helpers.py")]
    if ([p["path"] for p in stage["products"]] != [str(p) for p in expected_products]
            or [p["path"] for p in stage["modules"]] != [str(p) for p in expected_modules]):
        raise ValueError("Graph-stage output or scientific-module inventory differs")
    for pin in [*stage["products"], *stage["modules"]]:
        check(pin)


def run(path):
    plan_pin = record(path)
    plan = json.loads(path.read_text())
    source, science, output, repo = (Path(plan[k]) for k in ("source_root", "scientific_source", "output", "repo"))
    if (Path(__file__).resolve() != source / "benchmark_tools/continue_cap400_downstream.py"
            or (output / "results.json").exists()):
        raise ValueError("Wrong frozen runner or attempted namespace; no silent retry")
    if (os.environ.get("SLURM_CPUS_PER_TASK"), os.environ.get("SLURM_MEM_PER_NODE"),
            os.environ.get("SLURMD_NODENAME"), os.environ.get("CUDA_VISIBLE_DEVICES")) != (
            "32", "131072", "bizon", ""):
        raise ValueError("Require CPU-only 32-CPU 128-GiB bizon allocation")
    result = {"status": "validating_recovered_final_graph", "execution_plan": plan_pin,
        "job_id": os.environ.get("SLURM_JOB_ID"), "scope": plan["scope"], "scores": {},
        "defaults_promoted": False, "goal_completed": False, "host_before": host_activity(),
        "initial_hmm_search_rerun": False, "profile_expansion_rerun": False, "graph_clustering_rerun": False}
    write_json(output / "results.json", result)
    checked = [plan_pin, *plan["checked"], *plan["sources"]]
    try:
        for pin in checked:
            check(pin)
        parent, failure, admission, graph_plan = verify_graph(repo, plan["admission"])
        original_result = json.loads(Path(failure["parent_result"]["path"]).read_text())
        require_original_settings(plan, parent, admission, graph_plan, original_result)
        if runtime() != plan["runtime"] or plan["command"] != phylogeny_command(plan):
            raise ValueError("Changed original runtime or downstream command")
        env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(science),
            "PYTHONHASHSEED": "0", "CUDA_VISIBLE_DEVICES": "", "PYTHONFAULTHANDLER": "1",
            "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "32"}
        result["status"] = "original_refinement_and_candidates_running"
        write_json(output / "results.json", result)
        command = [sys.executable, str(source / "benchmark_tools/continue_cap400_downstream_worker.py"),
                   "--execution-plan", str(path)]
        result["graph_execution"] = timed_stage(command, output, "graph", science, env)
        if result["graph_execution"]["returncode"]:
            raise RuntimeError("Original refinement/candidate continuation failed; preserve attempt")
        stage_path = output / "graph_report.json"
        stage = json.loads(stage_path.read_text())
        validate_graph_stage(stage, plan, plan_pin)
        checked += [record(stage_path), *stage["products"], *stage["modules"]]
        result["graph"] = {"report": record(stage_path), "refined_groups": stage["refined_groups"],
                           "candidate_groups": stage["candidate_groups"]}
        result["status"] = "original_inferred_phylogeny_running"
        write_json(output / "results.json", result)
        result["execution"] = timed_stage(plan["command"], output, "native", science, env)
        if result["execution"]["returncode"]:
            raise RuntimeError("Original inferred phylogeny failed; preserve attempt")
        native_metrics = json.loads((output / "native_metrics.json").read_text())
        if (native_metrics.get("status") != "complete"
                or native_metrics.get("git", {}).get("commit") != REVISION
                or native_metrics.get("git", {}).get("dirty") is not False):
            raise ValueError("Native replay did not complete from original clean scientific revision")
        phylo = output / "native/orthohmm_phylogeny"
        prediction = output / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
        candidate = Path(argument(plan["command"], "--candidate-clusters"))
        result["coverage"] = validate_partition(candidate, prediction)
        result["outputs"] = [record(phylo / name) for name in PRODUCTS]
        native = replay_native(result, {"command": plan["command"]})
        if (native["genes"], native["species"]) != (984137, 78):
            raise ValueError("Recovered native inference changed the full input universe")
        checked += [record(output / "native_metrics.json"), *native["checked"]]
        result.update(status="native_verified_scoring_pending", native_summary=native["summary"],
                      native_checked=native["checked"], control=plan["control"])
        write_json(output / "results.json", result)
        result.update(score_qfo(repo, phylo, output, plan["scoring_id"]))
        scores, details, pins = readback_scores(repo, "QfO", result, native, output,
                                               "ohmm_hillclimb_r1_" + plan["scoring_id"])
        endpoint_vector("QfO", scores)
        result.update(scores=scores, details=details, scoring_readback=pins,
            deltas={k: value - plan["control"]["scores"][k] for k, value in scores.items()})
        checked += pins
        result["status"] = "original_cap400_cached_downstream_and_all_six_scores_verified"
    except Exception as error:
        result.update(status="failed_cap400_downstream_recovery", error_type=type(error).__name__, error=str(error))
    finally:
        try:
            for pin in checked:
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Runtime changed during downstream recovery")
            result["postflight_pins_verified"] = True
        except Exception as error:
            result.update(status="failed_cap400_downstream_recovery", postflight_error=str(error))
        result["host_after"] = host_activity()
        write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    for name in ("repo", "source", "output"):
        prep.add_argument("--" + name, required=True, type=Path)
    worker = sub.add_parser("run")
    worker.add_argument("--execution-plan", required=True, type=Path)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve())
              if args.command == "prepare" else run(args.execution_plan.resolve()))
    print(result["status"])
    raise SystemExit(0 if result["status"] in {"prepared_unrun_cap400_downstream",
        "original_cap400_cached_downstream_and_all_six_scores_verified"} else 1)

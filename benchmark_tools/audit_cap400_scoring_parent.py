"""Authenticate completed cap400 native inference after a scoring import failure."""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.continue_cap400_downstream import validate_graph_stage
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.report_cap400_downstream import load_plan
from benchmark_tools.report_hillclimb import replay_native
from benchmark_tools.run_membership_hillclimb import argument, write_json


def require_failure(report, rows):
    expected = {"JobIDRaw": "24350", "State": "FAILED", "ExitCode": "1:0",
                "AllocCPUS": "32", "NodeList": "bizon"}
    if rows != [expected] or report.get("job_id") != "24350":
        raise ValueError("Require the unique terminal original downstream failure")
    if (report.get("status") != "failed_cap400_downstream_recovery"
            or report.get("error_type") != "ModuleNotFoundError"
            or report.get("error") != "No module named 'qfo_benchmark'"
            or report.get("scores") != {} or report.get("postflight_pins_verified") is not True
            or any(report.get(key) is not False for key in (
                "initial_hmm_search_rerun", "profile_expansion_rerun", "graph_clustering_rerun"))
            or "postflight_error" in report):
        raise ValueError("Not the specific post-native scoring-import failure")
    return expected


def require_native_execution(root, plan, report, metrics):
    graph = [sys.executable, str(Path(plan["source_root"]) / "benchmark_tools/continue_cap400_downstream_worker.py"),
             "--execution-plan", str(root / "execution_plan.json")]
    if (report["graph_execution"].get("command") != graph
            or report["graph_execution"].get("returncode") != 0
            or report["execution"].get("command") != plan["command"]
            or report["execution"].get("returncode") != 0):
        raise ValueError("Graph or native inference did not finish with the original command")
    parameters = {"species_tree_mode": "infer", "species_tree": None,
        "aligner": plan["aligner"], "tree_builder": plan["tree_builder"],
        "root_duplication_rule": "species_overlap", "pair_orthology_rule": "positive_paralogy",
        "species_tree_rooting": "min_variance", "cpu": 32,
        "satellite_membership_policy": "high_confidence_pair",
        "explicit_unconstrained_ablation": False, "checkpoint_source": None}
    source = Path(plan["scientific_source"]) / "benchmark_tools/replay_phylogeny.py"
    if (metrics.get("status") != "complete" or metrics.get("source") != record(source)
            or metrics.get("git", {}).get("commit") != plan["scientific_revision"]
            or metrics.get("git", {}).get("dirty") is not False
            or metrics.get("parameters") != parameters):
        raise ValueError("Native metrics do not prove unchanged completed original inference")
    expected = {"candidate_clusters": record(argument(plan["command"], "--candidate-clusters")),
        "membership_constraints": record(argument(plan["command"], "--membership-constraints")),
        "fasta_directory": plan["fasta"]}
    if any(metrics.get("input", {}).get(key) != value for key, value in expected.items()):
        raise ValueError("Native inference used different candidate, constraints or FASTAs")


def audit(repo, root):
    path = root / "results.json"
    producer = record(path)
    report = json.loads(path.read_text())
    accounting = subprocess.check_output(["sacct", "-j", "24350", "-X", "-P", "-o",
        "JobIDRaw,State,ExitCode,AllocCPUS,NodeList"], text=True)
    scheduler = require_failure(report, list(csv.DictReader(io.StringIO(accounting), delimiter="|")))
    plan, pins, failed_fresh = load_plan(repo, root, report)
    metrics_path = root / "native_metrics.json"
    metrics = json.loads(metrics_path.read_text())
    require_native_execution(root, plan, report, metrics)
    stage_pin = report["graph"]["report"]
    if Path(stage_pin["path"]).resolve() != root / "graph_report.json":
        raise ValueError("Graph stage belongs to another execution")
    check(stage_pin)
    stage = json.loads(Path(stage_pin["path"]).read_text())
    validate_graph_stage(stage, plan, report["execution_plan"])
    if any(report["graph"][key] != stage[key] for key in ("refined_groups", "candidate_groups")):
        raise ValueError("Reported graph-stage counts differ")
    native = replay_native(report, {"command": plan["command"]})
    if ((native["genes"], native["species"]) != (984137, 78)
            or native["summary"] != report.get("native_summary")
            or native["checked"] != report.get("native_checked")):
        raise ValueError("Retained native content or full universe differs from parent readback")
    work = repo / "qfo_benchmark/w" / ("hc1_" + plan["scoring_id"])
    if work.exists() or (root / "qfo_scoring").exists():
        raise ValueError("Parent unexpectedly contains scoring work")
    checked = [producer, *pins, stage_pin, record(metrics_path), *stage["products"],
               *stage["modules"], *native["checked"], *plan["runtime"]["files"]]
    checked = list({pin["path"]: pin for pin in checked}.values())
    for pin in checked:
        check(pin)
    return {"status": "completed_original_cap400_native_reverified_scoring_only_pending",
        "producer": producer, "execution_plan": report["execution_plan"], "scheduler": scheduler,
        "scientific_revision": plan["scientific_revision"], "source": record(__file__),
        "native": native, "checked": checked, "original_failed_fresh_attempt": failed_fresh,
        "execution": report["execution"], "scores": {}, "native_inference_rerun": False,
        "scope": "Post-native scoring-import failure authenticated. Original products untouched; no accuracy admission or fresh full-pipeline timing.",
        "defaults_promoted": False, "goal_completed": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "parent-root", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.repo.resolve(), args.parent_root.resolve())
    write_json(args.output.resolve(), result)
    print(result["status"])

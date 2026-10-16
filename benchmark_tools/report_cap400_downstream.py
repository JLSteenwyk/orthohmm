"""Admit recovered cap400 scores while preserving the failed fresh attempt."""

import csv
import io
import json
from pathlib import Path
import subprocess
import sys

from benchmark_tools.continue_cap400_downstream import (
    require_original_settings, validate_graph_stage, verify_graph,
)
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.recover_full_qfo_search import short_id
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native, snapshot
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import argument

SUCCESS = "original_cap400_cached_downstream_and_all_six_scores_verified"


def terminal(report):
    if report.get("job_id") != "24350":
        raise ValueError("Wrong original cap400 downstream producer job")
    text = subprocess.check_output(["sacct", "-j", "24350", "-X", "-P", "-o",
        "JobIDRaw,State,ExitCode,AllocCPUS,NodeList"], text=True)
    rows = list(csv.DictReader(io.StringIO(text), delimiter="|"))
    expected = {"JobIDRaw": "24350", "State": "COMPLETED", "ExitCode": "0:0",
                "AllocCPUS": "32", "NodeList": "bizon"}
    if rows != [expected]:
        raise ValueError("Cap400 downstream producer is not terminal successful on planned resources")
    return expected


def capture(path, directory):
    pin = record(path)
    report = json.loads(path.read_text())
    check(pin)
    if report.get("status") != SUCCESS:
        return snapshot(path, directory, "cap400_recovery.json")
    scheduler = terminal(report)
    reference = directory / "cap400_recovery_reference.json"
    from benchmark_tools.run_membership_hillclimb import write_json
    write_json(reference, {"producer": pin, "scheduler": scheduler,
        "scope": "Immutable terminal producer; complete verification inventories remain in its pinned report and execution plan."})
    return report, {"snapshot": pin, "reference": record(reference)}


def load_plan(repo, root, report):
    check(report["execution_plan"])
    if Path(report["execution_plan"]["path"]).resolve() != root / "execution_plan.json":
        raise ValueError("Cap400 producer plan belongs to another namespace")
    plan = json.loads(Path(report["execution_plan"]["path"]).read_text())
    if Path(plan["repo"]).resolve() != repo or Path(plan["output"]).resolve() != root:
        raise ValueError("Cap400 recovery belongs to another repository or output")
    pins = [report["execution_plan"], *plan["checked"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in pins:
        check(pin)
    if runtime() != plan["runtime"]:
        raise ValueError("Original scientific runtime changed")
    parent, failure, admission, graph_plan = verify_graph(repo, plan["admission"])
    original = json.loads(Path(failure["parent_result"]["path"]).read_text())
    require_original_settings(plan, parent, admission, graph_plan, original)
    failed = {"dataset": "QfO", "status": "failed", "scores": {},
        "producer_status": original["status"], "producer_job_id": original["job_id"],
        "source_snapshot": failure["parent_result"], "error": original["error"],
        "error_type": original["error_type"], "resources": original["execution"],
        "execution_kind": "failed_fresh_original_arrival_cap400_final_clustering",
        "scope": "Original failed fresh run retained; no scores from its incomplete partition."}
    return plan, pins, failed


def require_completed_execution(root, plan, report):
    if (report.get("status") != SUCCESS or report.get("postflight_pins_verified") is not True
            or any(report.get(key) is not False for key in (
                "initial_hmm_search_rerun", "profile_expansion_rerun", "graph_clustering_rerun"))):
        raise ValueError("Require completed verified original cap400 downstream recovery")
    expected = [sys.executable, str(Path(plan["source_root"]) / "benchmark_tools/continue_cap400_downstream_worker.py"),
                "--execution-plan", str(root / "execution_plan.json")]
    if (report["graph_execution"].get("command") != expected
            or report["graph_execution"].get("returncode") != 0
            or report["execution"].get("command") != plan["command"]
            or report["execution"].get("returncode") != 0
            or plan["scoring_id"] != short_id(Path(plan["repo"]), root, 1)):
        raise ValueError("Actual cap400 execution or scoring namespace differs")


def control_readback(repo, plan, report, original):
    claimed = original["control"]
    if plan["control"] != claimed or report["control"] != claimed:
        raise ValueError("Recovery changed the matched original control")
    pin = claimed["report"]
    check(pin)
    root = Path(pin["path"]).parent
    control = json.loads(Path(pin["path"]).read_text())
    if control.get("status") != "cached_native_control_and_all_six_scores_verified":
        raise ValueError("Original control is not independently verified")
    check(control["execution_plan"])
    control_plan = json.loads(Path(control["execution_plan"]["path"]).read_text())
    native = replay_native(control, {"command": control_plan["command"]})
    preflight = root / "scoring_preflight.json"
    participant = json.loads(preflight.read_text())["stage"]["participant"]
    scores, details, pins = readback_scores(repo, "QfO", control, native, root, participant)
    if scores != claimed["scores"] or (native["genes"], native["species"]) != (984137, 78):
        raise ValueError("Matched original control content or raw scores differ")
    return {"report": pin, "scores": endpoint_vector("QfO", scores), "details": details,
        "products": native["products"], "coverage": {k: native[k] for k in ("genes", "species", "groups", "pairs")},
        "source_revision": plan["scientific_revision"]}, [record(preflight), *native["checked"], *pins]


def cell(repo, root, plan, report):
    require_completed_execution(root, plan, report)
    terminal(report)
    parent, failure, admission, graph_plan = verify_graph(repo, plan["admission"])
    original = json.loads(Path(failure["parent_result"]["path"]).read_text())
    require_original_settings(plan, parent, admission, graph_plan, original)
    if Path(report["graph"]["report"]["path"]).resolve() != root / "graph_report.json":
        raise ValueError("Graph-stage report belongs to another recovery")
    check(report["graph"]["report"])
    stage = json.loads(Path(report["graph"]["report"]["path"]).read_text())
    validate_graph_stage(stage, plan, report["execution_plan"])
    if any(report["graph"][k] != stage[k] for k in ("refined_groups", "candidate_groups")):
        raise ValueError("Recovered graph-stage counts differ")
    metrics_path = root / "native_metrics.json"
    metrics = json.loads(metrics_path.read_text())
    source = Path(plan["scientific_source"]) / "benchmark_tools/replay_phylogeny.py"
    parameters = {"species_tree_mode": "infer", "species_tree": None,
        "aligner": plan["aligner"], "tree_builder": plan["tree_builder"],
        "root_duplication_rule": "species_overlap", "pair_orthology_rule": "positive_paralogy",
        "species_tree_rooting": "min_variance", "cpu": 32,
        "satellite_membership_policy": "high_confidence_pair",
        "explicit_unconstrained_ablation": False, "checkpoint_source": None}
    if (metrics.get("status") != "complete" or metrics.get("source") != record(source)
            or metrics.get("git", {}).get("commit") != plan["scientific_revision"]
            or metrics.get("git", {}).get("dirty") is not False
            or metrics.get("parameters") != parameters):
        raise ValueError("Recovered native metrics have wrong source, scope or resources")
    if (metrics.get("input", {}).get("candidate_clusters") != record(argument(plan["command"], "--candidate-clusters"))
            or metrics["input"].get("membership_constraints") != record(argument(plan["command"], "--membership-constraints"))
            or metrics["input"].get("fasta_directory") != plan["fasta"]):
        raise ValueError("Native candidate/constraint/FASTA inputs differ from recovery plan")
    native = replay_native(report, {"command": plan["command"]})
    if (native["genes"], native["species"]) != (984137, 78):
        raise ValueError("Recovered cap400 native gene universe differs")
    scores, details, pins = readback_scores(repo, "QfO", report, native, root,
                                           "ohmm_hillclimb_r1_" + plan["scoring_id"])
    matched, control_pins = control_readback(repo, plan, report, original)
    for key, value in scores.items():
        if report["deltas"][key] != value - original["control"]["scores"][key]:
            raise ValueError("Reported cap400 delta differs from original matched control")
    checked = [*plan["checked"], *plan["sources"], *native["checked"], *pins, *control_pins,
        record(metrics_path), report["graph"]["report"], *stage["products"], *stage["modules"]]
    for pin in checked:
        check(pin)
    graph_report = json.loads(Path(admission["producer"]["path"]).read_text())
    compact = [report["execution_plan"], record(root / "results.json"), plan["admission"],
        admission["execution_plan"], admission["producer"], failure["parent_result"],
        record(metrics_path), report["graph"]["report"], *stage["products"], *stage["modules"],
        *native["checked"], *report["metric_files"], matched["report"], *original["control"]["checked"]]
    return {"dataset": "QfO", "status": "independent_content_readback_verified",
        "scores": endpoint_vector("QfO", scores), "details": details, "products": native["products"],
        "checked": list({p["path"]: p for p in compact}.values()),
        "verification_record_count": len({p["path"] for p in checked}),
        "verification_scope": "Full producer and raw-score inventories rehashed; complete inventories remain in pinned plans/reports rather than duplicated in this table.",
        "coverage": {k: native[k] for k in ("genes", "species", "groups", "pairs")},
        "source_revision": plan["scientific_revision"], "matched_control": matched,
        "original_failed_attempt": failure["parent_result"],
        "resources": {"failed_fresh_execution": original["execution"],
            "cached_final_graph_recovery": graph_report["execution"],
            "cached_downstream_graph": report["graph_execution"], "downstream_native": report["execution"]},
        "execution_kind": "cached_final_graph_recovery_of_original_arrival_cap400",
        "scope": "Original failed fresh run preserved. Cached graph/downstream recovery is not uninterrupted full-pipeline timing. FAS remains unseeded."}

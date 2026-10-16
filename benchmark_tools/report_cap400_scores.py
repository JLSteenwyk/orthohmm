"""Admit scoring-only cap400 recovery without rewriting either failed attempt."""

import csv
import io
import json
from pathlib import Path
import subprocess

from benchmark_tools.audit_cap400_scoring_parent import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.recover_cap400_scores import (
    SUCCESS, read_admission, require_admission, verify_helpers,
)
from benchmark_tools.recover_full_qfo_search import short_id
from benchmark_tools.report_cap400_downstream import control_readback
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, snapshot
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.summarize_matched_resources import verbose_time

SUBMISSION = "benchmark_tools/results/cap400_scoring_only_submission_20261016.json"


def terminal(report):
    if report.get("job_id") != "24355":
        raise ValueError("Wrong cap400 scoring-only producer job")
    text = subprocess.check_output(["sacct", "-j", "24355", "-X", "-P", "-o",
        "JobIDRaw,State,ExitCode,AllocCPUS,NodeList"], text=True)
    rows = list(csv.DictReader(io.StringIO(text), delimiter="|"))
    expected = {"JobIDRaw": "24355", "State": "COMPLETED", "ExitCode": "0:0",
                "AllocCPUS": "8", "NodeList": "bizon"}
    if rows != [expected]:
        raise ValueError("Scoring-only producer is not uniquely terminal successful on planned resources")
    return expected


def capture(path, directory):
    pin = record(path)
    report = json.loads(path.read_text())
    check(pin)
    if report.get("status") != SUCCESS:
        return snapshot(path, directory, "cap400_scoring.json")
    scheduler = terminal(report)
    reference = directory / "cap400_scoring_reference.json"
    write_json(reference, {"producer": pin, "scheduler": scheduler,
        "scope": "Immutable terminal scoring-only producer; full inventories remain in pinned plans and native audit."})
    return report, {"snapshot": pin, "reference": record(reference)}


def load_plan(repo, root, report):
    launch_pin = record(repo / SUBMISSION)
    launch = json.loads(Path(launch_pin["path"]).read_text())
    if (report.get("execution_plan") != launch["execution_plan"]
            or report.get("job_id") != launch["job_id"]
            or Path(report["execution_plan"]["path"]).resolve() != root / "execution_plan.json"):
        raise ValueError("Scoring-only report differs from authenticated launch")
    check(report["execution_plan"])
    plan = json.loads(Path(report["execution_plan"]["path"]).read_text())
    if (Path(plan["repo"]).resolve() != repo or Path(plan["output"]).resolve() != root
            or plan["source_root"] != launch["source_root"]
            or plan["source_revision"] != launch["source_revision"]
            or report.get("source_revision") != plan["source_revision"]
            or report.get("scientific_revision") != plan["scientific_revision"]
            or plan["scientific_revision"] != launch["scientific_revision"]
            or plan.get("native_inference_rerun") is not False
            or short_id(repo, root, 1) != plan["scoring_id"]):
        raise ValueError("Recovery plan changed executor, original science or namespace")
    pins = [launch_pin, report["execution_plan"], plan["admission"], plan["parent_plan"],
            plan["parent_report"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in pins:
        check(pin)
    admission, evidence = read_admission(plan["admission"])
    parent = json.loads(Path(plan["parent_plan"]["path"]).read_text())
    if (evidence["execution_plan"] != plan["parent_plan"] or evidence["producer"] != plan["parent_report"]
            or parent["output"] != plan["parent_root"] or parent["repo"] != str(repo)
            or runtime() != plan["runtime"] or parent["runtime"] != plan["runtime"]
            or verify_helpers(Path(plan["source_root"]), parent, plan["source_revision"]) != plan["sources"]):
        raise ValueError("Recovery provenance or original scientific runtime changed")
    failed_parent = json.loads(Path(plan["parent_report"]["path"]).read_text())
    failed = {"dataset": "QfO", "status": "failed", "scores": {},
        "producer_status": failed_parent["status"], "producer_job_id": "24350",
        "source_snapshot": plan["parent_report"], "error": failed_parent["error"],
        "error_type": failed_parent["error_type"],
        "resources": {"graph": failed_parent["graph_execution"], "native": failed_parent["execution"]},
        "execution_kind": "failed_scoring_import_after_completed_original_native",
        "scope": "Terminal failed job retained; native stage succeeded but supplied no scores."}
    pins += [admission["audit"], *evidence["checked"]]
    return plan, pins, [evidence["original_failed_fresh_attempt"], failed]


def require_completed(report):
    if (report.get("status") != SUCCESS or report.get("postflight_pins_verified") is not True
            or report.get("native_inference_rerun") is not False
            or report.get("defaults_promoted") is not False or report.get("goal_completed") is not False
            or "error" in report or "postflight_error" in report):
        raise ValueError("Require successfully verified scoring-only recovery")


def cell(repo, root, plan, report):
    require_completed(report)
    scheduler = terminal(report)
    admission, evidence = read_admission(plan["admission"])
    fresh = audit(repo, Path(plan["parent_root"]))
    require_admission(admission, fresh)
    native = fresh["native"]
    if (native != evidence["native"] or report.get("native_products") != native["products"]
            or report.get("native_summary") != native["summary"]
            or report.get("coverage") != {k: native[k] for k in ("genes", "species", "groups", "pairs")}
            or report.get("native_audit") != admission["audit"]
            or report.get("converter") != record(Path(plan["source_root"]) / "qfo_benchmark/og_to_pairwise.py")):
        raise ValueError("Recovered scores do not belong to the authenticated original native products")
    scores, details, score_pins = readback_scores(repo, "QfO", report["scoring"], native, root,
                                               "ohmm_hillclimb_r1_" + plan["scoring_id"])
    parent = json.loads(Path(plan["parent_plan"]["path"]).read_text())
    failed_parent = json.loads(Path(plan["parent_report"]["path"]).read_text())
    original_pin = evidence["original_failed_fresh_attempt"]["source_snapshot"]
    original = json.loads(Path(original_pin["path"]).read_text())
    matched, control_pins = control_readback(repo, parent, failed_parent, original)
    if (scores != report["scores"] or report.get("deltas") != {
            key: value - original["control"]["scores"][key] for key, value in scores.items()}):
        raise ValueError("Recovered cap400 scores or matched-control deltas differ")
    time_pin = record(root / "recovery.time.txt")
    resources = verbose_time(Path(time_pin["path"]).read_text())
    graph_admission = json.loads(Path(parent["admission"]["path"]).read_text())
    graph_report = json.loads(Path(graph_admission["producer"]["path"]).read_text())
    checked = [plan["admission"], admission["audit"], report["execution_plan"],
        record(root / "results.json"), plan["parent_plan"], plan["parent_report"], original_pin,
        time_pin, report["converter"], graph_admission["producer"],
        *fresh["checked"], *score_pins, *control_pins]
    for pin in checked:
        check(pin)
    compact = [plan["admission"], admission["audit"], report["execution_plan"], record(root / "results.json"),
        plan["parent_plan"], plan["parent_report"], original_pin, time_pin, report["converter"],
        graph_admission["producer"], matched["report"], *native["checked"], *score_pins]
    return {"dataset": "QfO", "status": "independent_content_readback_verified",
        "scores": endpoint_vector("QfO", scores), "details": details, "products": native["products"],
        "checked": list({p["path"]: p for p in compact}.values()),
        "verification_record_count": len({p["path"] for p in checked}),
        "coverage": {k: native[k] for k in ("genes", "species", "groups", "pairs")},
        "source_revision": plan["scientific_revision"], "matched_control": matched,
        "resources": {"failed_fresh_execution": original["execution"],
            "cached_final_graph_recovery": graph_report["execution"],
            "cached_downstream_graph": failed_parent["graph_execution"],
            "completed_native_stage_with_failed_scoring_import": failed_parent["execution"],
            "scoring_only_recovery": {"scheduler": scheduler, "gnu_time": resources, "source": time_pin}},
        "execution_kind": "scoring_only_recovery_of_original_arrival_cap400_native",
        "scope": "Both failed attempts retained. Original native products and all six raw scores independently read back. Separate recovery costs are not uninterrupted full-pipeline timing. FAS remains unseeded."}

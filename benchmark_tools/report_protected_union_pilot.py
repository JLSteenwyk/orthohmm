"""Independently admit the fresh OrthoBench protected-union pilot, not partial scores."""

import csv
import io
import json
import math
from pathlib import Path
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, snapshot
from benchmark_tools.run_full_search_hillclimb import runtime, verify_input_inventory
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_protected_union_pilot import (
    PROTOCOL, SLURM, SUCCESS, committed, pilot_gate, prerequisites, scientific_source,
    validate_plan, validate_variant,
)
from benchmark_tools.run_three_kingdoms_membership_supplement import sources
from benchmark_tools.summarize_matched_resources import verbose_time

SUBMISSION = "benchmark_tools/results/protected_union_pilot_submission_20261016.json"
LAUNCH_REVISION = "c7a9c36a82a6983936e7c3747a1eedcd31447d03"
REVISION = "d650d6853869238a738b31367688018158d49625"
KEY = "protected_union:protected_union"
LABEL = "OrthoHMM protected candidate union (experimental)"


def terminal(report):
    if report.get("job_id") != "24363":
        raise ValueError("Wrong protected-union OrthoBench pilot job")
    text = subprocess.check_output(["sacct", "-j", "24363", "-X", "-P", "-o",
        "JobIDRaw,State,ExitCode,AllocCPUS,NodeList"], text=True)
    rows = list(csv.DictReader(io.StringIO(text), delimiter="|"))
    expected = {"JobIDRaw": "24363", "State": "COMPLETED", "ExitCode": "0:0",
                "AllocCPUS": "16", "NodeList": "bizon"}
    if rows != [expected]:
        raise ValueError("Pilot is not uniquely terminal successful on planned resources")
    return expected


def capture(path, directory):
    pin = record(path)
    report = json.loads(path.read_text())
    check(pin)
    if report.get("status") != SUCCESS:
        return snapshot(path, directory, "protected_pilot.json")
    reference = directory / "protected_pilot_reference.json"
    write_json(reference, {"producer": pin, "scheduler": terminal(report),
        "scope": "Terminal fresh full-proteome OrthoBench pilot; gate result does not imply default promotion."})
    return report, {"snapshot": pin, "reference": record(reference)}


def load_plan(repo, root, report):
    launch, launch_pin = committed(repo, repo, LAUNCH_REVISION, SUBMISSION)
    if (launch.get("job_id") != "24363" or launch.get("source_revision") != REVISION
            or report.get("execution_plan") != launch["execution_plan"]
            or report.get("job_id") != launch["job_id"]
            or report.get("source_revision") != REVISION
            or Path(report["execution_plan"]["path"]).resolve() != root / "execution_plan.json"):
        raise ValueError("Pilot report differs from committed launch")
    check(report["execution_plan"])
    plan = json.loads(Path(report["execution_plan"]["path"]).read_text())
    source = Path(plan["source_root"])
    if (Path(plan["repo"]).resolve() != repo or Path(plan["output"]).resolve() != root
            or plan["source_root"] != launch["source_root"] or plan["source_revision"] != REVISION
            or plan.get("defaults_promoted") is not False or plan.get("goal_completed") is not False
            or subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip() != REVISION):
        raise ValueError("Pilot source, namespace or promotion state differs from launch")
    pins = [launch_pin, report["execution_plan"], launch["working_tests"], launch["frozen_tests"],
            *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in pins:
        check(pin)
    expected_sources = [*sources(source), record(source / "benchmark_tools" / PROTOCOL),
                        record(source / "benchmark_tools" / SLURM)]
    if expected_sources != plan["sources"] or runtime() != plan["runtime"]:
        raise ValueError("Frozen pilot source inventory or scientific runtime changed")
    parent, fresh, control, parent_pins, accounting = prerequisites(repo, source, REVISION)
    if control != plan["control"] or accounting != plan["dependency_accounting"]:
        raise ValueError("Pilot matched control or terminal prerequisite accounting changed")
    validate_plan(plan, parent)
    pins += [*parent_pins, *scientific_source(source, fresh)]
    verify_input_inventory(plan["case"])
    return plan, list({p["path"]: p for p in pins}.values()), []


def require_completed(report, plan):
    if (report.get("status") != SUCCESS or report.get("postflight_pins_verified") is not True
            or report.get("defaults_promoted") is not False or report.get("goal_completed") is not False
            or report.get("source_revision") != plan["source_revision"]
            or "error" in report or "postflight_error" in report):
        raise ValueError("Require successfully verified full OrthoBench pilot")


def cell(repo, root, plan, report):
    require_completed(report, plan)
    scheduler = terminal(report)
    native, routes, native_pins = validate_variant(plan)
    coverage = {k: native[k] for k in ("genes", "species", "groups", "pairs")}
    if ((native["genes"], native["species"]) != (251378, 12)
            or report.get("coverage") != coverage or report.get("native_products") != native["products"]
            or report.get("native_summary") != native["summary"] or report.get("route_evidence") != routes
            or report.get("control") != plan["control"]):
        raise ValueError("Pilot report differs from complete native outputs, HMM routes or matched control")
    execution = report["execution"]
    folder = Path(plan["case"]["output"])
    time_pin, log_pin = record(folder / "native.time.txt"), record(folder / "native.log")
    resources = verbose_time(Path(time_pin["path"]).read_text())
    if (execution.get("returncode") != 0 or execution.get("command") != plan["case"]["command"]
            or execution.get("time") != time_pin or execution.get("log") != log_pin
            or execution.get("resources") != resources):
        raise ValueError("Pilot native command, logs or resource record differs")
    wall = report.get("total_pilot_wall_s")
    if isinstance(wall, bool) or not isinstance(wall, (int, float)) or not math.isfinite(wall) or wall <= 0:
        raise ValueError("Require finite positive pilot wall time")
    scores, details, score_pins = readback_scores(repo, "OrthoBench", report["scoring"], native,
                                               folder, "unused_for_orthobench")
    gate = pilot_gate(scores, details, plan["control"])
    if (scores != report.get("scores") or details != report.get("details") or gate != report.get("gate")
            or report["scoring"].get("scores") != scores
            or report["scoring"].get("details") != details["family_details"]):
        raise ValueError("Independent pilot scores, precision/recall or precommitted gate differs")
    checked = [report["execution_plan"], record(root / "results.json"), time_pin, log_pin,
               *native_pins, *score_pins, *plan["control"]["checked"]]
    checked = list({p["path"]: p for p in checked}.values())
    for pin in checked:
        check(pin)
    return {"dataset": "OrthoBench", "status": "independent_content_readback_verified",
        "scores": endpoint_vector("OrthoBench", scores), "details": details,
        "coverage": coverage, "products": native["products"], "checked": checked,
        "source_revision": REVISION, "gate": gate, "route_evidence": routes,
        "matched_control": plan["control"],
        "resources": {"scheduler": scheduler, "native": execution, "total_pilot_wall_s": wall,
            "host_before": report["host_before"], "host_after": report["host_after"]},
        "execution_kind": "fresh_full_proteome_protected_union_with_inferred_phylogeny",
        "scope": "Fresh OrthoBench inference and independent score readback. Cached matched control has different timing scope; shared-host timings are descriptive. QfO/Three Kingdoms are absent until separately prespecified, executed and verified. A passing pilot gate is not an eight-endpoint improvement, equivalence or default promotion."}


def require_existing_control(methods, admitted):
    controls = [m for m in methods if m["key"] == "canonical_collection:canonical"]
    expected = admitted["matched_control"]
    if len(controls) != 1:
        raise ValueError("Existing canonical control is missing or ambiguous")
    cells = [c for c in controls[0]["cells"] if c["dataset"] == "OrthoBench"
             and c["status"] == "independent_content_readback_verified"]
    if (len(cells) != 1 or {k: v for k, v in controls[0]["scores"].items() if k.startswith("OrthoBench:")}
            != endpoint_vector("OrthoBench", expected["scores"])
            or any(cells[0][k] != expected[k] for k in ("details", "coverage", "products", "source_revision"))):
        raise ValueError("Existing canonical control differs from independently reverified matched control")

"""Admit a proven lossless shared-QfO export view, retaining the failed producer."""

import csv
import io
import json
from pathlib import Path
import subprocess

from benchmark_tools.continue_shared_qfo import REVISION, original_skip
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.recover_cap400_scores import verify_helpers
from benchmark_tools.recover_full_qfo_search import short_id
from benchmark_tools.recover_shared_qfo_export import (
    DISTINCT, DUPLICATE, DUPLICATE_LINE, PROTOCOL, ROWS, SLURM, SUCCESS, parent_context,
)
from benchmark_tools.report_additional_hillclimb import shared_native
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, snapshot
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json
from benchmark_tools.shared_qfo_export import verify_projection
from benchmark_tools.summarize_matched_resources import verbose_time

SUBMISSION = "benchmark_tools/results/shared_qfo_export_recovery_submission_20261016.json"


def terminal(report):
    if report.get("job_id") != "24360":
        raise ValueError("Wrong shared-QfO export-recovery job")
    text = subprocess.check_output(["sacct", "-j", "24360", "-X", "-P", "-o",
        "JobIDRaw,State,ExitCode,AllocCPUS,NodeList"], text=True)
    rows = list(csv.DictReader(io.StringIO(text), delimiter="|"))
    expected = {"JobIDRaw": "24360", "State": "COMPLETED", "ExitCode": "0:0",
                "AllocCPUS": "8", "NodeList": "bizon"}
    if rows != [expected]:
        raise ValueError("Export recovery is not uniquely terminal successful on planned resources")
    return expected


def capture(path, directory):
    pin = record(path)
    report = json.loads(path.read_text())
    check(pin)
    if report.get("status") != SUCCESS:
        return snapshot(path, directory, "shared_export.json")
    reference = directory / "shared_export_reference.json"
    write_json(reference, {"producer": pin, "scheduler": terminal(report),
        "scope": "Terminal derived-export recovery; original failed producer remains independently pinned."})
    return report, {"snapshot": pin, "reference": record(reference)}


def load_plan(repo, root, report):
    launch_pin = record(repo / SUBMISSION)
    launch = json.loads(Path(launch_pin["path"]).read_text())
    if (report.get("execution_plan") != launch["execution_plan"]
            or report.get("job_id") != launch["job_id"]
            or Path(report["execution_plan"]["path"]).resolve() != root / "execution_plan.json"):
        raise ValueError("Export-recovery report differs from authenticated launch")
    check(report["execution_plan"])
    plan = json.loads(Path(report["execution_plan"]["path"]).read_text())
    if (Path(plan["repo"]).resolve() != repo or Path(plan["output"]).resolve() != root
            or plan["source_root"] != launch["source_root"]
            or plan["source_revision"] != launch["source_revision"]
            or report.get("source_revision") != plan["source_revision"]
            or report.get("scientific_revision") != plan["scientific_revision"]
            or plan["scientific_revision"] != launch["scientific_revision"]
            or plan["scientific_revision"] != REVISION
            or plan.get("native_inference_rerun") is not False
            or short_id(repo, root, 1) != plan["scoring_id"]):
        raise ValueError("Export recovery changed executor, original science or namespace")
    pins = [launch_pin, report["execution_plan"], plan["failure"],
            *plan["sources"], *plan["runtime"]["files"]]
    for pin in pins:
        check(pin)
    original, parent, failed, parent_pins = parent_context(repo, plan["failure"])
    if (runtime() != plan["runtime"] or original["runtime"] != plan["runtime"]
            or verify_helpers(Path(plan["source_root"]), original, plan["source_revision"],
                              protocol=PROTOCOL, slurm=SLURM) != plan["sources"]):
        raise ValueError("Export-recovery provenance or original scientific runtime changed")
    skip_pin = failed["original_skip"]
    check(skip_pin)
    if Path(skip_pin["path"]).resolve() != Path(original["output"]) / "original_skip_snapshot.json":
        raise ValueError("Original skipped attempt belongs to another continuation")
    archived = json.loads(Path(skip_pin["path"]).read_text())
    skipped = original_skip(archived, parent, original["parent_plan"])
    failures = [{**skipped, "scores": {}, "source_snapshot": skip_pin,
                 "producer_job_id": archived["job_id"]},
        {"dataset": "QfO", "status": "failed", "scores": {},
         "producer_status": failed["status"], "producer_job_id": "24247",
         "source_snapshot": record(Path(original["output"]) / "results.json"),
         "error": failed["error"], "error_type": failed["error_type"],
         "resources": failed["execution"],
         "execution_kind": "failed_native_export_after_completed_original_inference",
         "scope": "Original native command exited zero but strict validation rejected a redundant pair; no original scores admitted."}]
    return plan, [*pins, *parent_pins, skip_pin], failures


def require_completed(report):
    if (report.get("status") != SUCCESS or report.get("postflight_pins_verified") is not True
            or report.get("native_inference_rerun") is not False
            or report.get("defaults_promoted") is not False or report.get("goal_completed") is not False
            or "error" in report or "postflight_error" in report):
        raise ValueError("Require successfully verified lossless export recovery")


def export_readback(root, report, original, parent):
    manifest_pin = report["normalization"]
    check(manifest_pin)
    if Path(manifest_pin["path"]).resolve() != root / "normalization.json":
        raise ValueError("Normalization manifest belongs to another recovery")
    manifest = json.loads(Path(manifest_pin["path"]).read_text())
    proof = manifest["normalization"]
    view = root / "canonical_export"
    folder = Path(original["case"]["output"])
    phylo = folder / "output/orthohmm_phylogeny"
    if (manifest.get("status") != "explicit_lossless_canonical_export_view_verified"
            or manifest.get("native_inference_rerun") is not False
            or proof.get("duplicate") != list(DUPLICATE)
            or proof.get("duplicate_confidence") != "medium"
            or proof.get("original_rows") != ROWS or proof.get("distinct_rows") != DISTINCT
            or any(proof["skips"][name]["line"] != DUPLICATE_LINE for name in ("pairs", "confidence"))
            or proof["original"] != {name: record(phylo / product)
                                      for name, product in zip(("pairs", "confidence"), PRODUCTS[1:])}
            or proof["derived"] != {name: record(view / product)
                                     for name, product in zip(("pairs", "confidence"), PRODUCTS[1:])}
            or manifest["original_completed_native_metadata"] != record(folder / "metrics.json")):
        raise ValueError("Normalization differs from the exact planned original export")
    projection = verify_projection(proof)
    if projection != manifest["independent_projection"]:
        raise ValueError("Independent byte projection differs from recorded proof")
    native = shared_native(original["case"], parent, pair_export=proof)
    if (native["status"] != "derived_native_content_verified"
            or (native["genes"], native["species"], native["groups"], native["pairs"])
            != (984137, 78, 370016, DISTINCT)):
        raise ValueError("Derived view differs from complete original native coverage")
    copied = []
    for name in (PRODUCTS[0], "provenance_manifest.json", "reconciliation_summary.json", "species_tree.rooted.nwk"):
        before, after = record(phylo / name), record(view / name)
        if (before["bytes"], before["sha256"]) != (after["bytes"], after["sha256"]):
            raise ValueError("Derived metadata or root groups differ from original")
        copied += [before, after]
    if copied != manifest["unchanged_copied_products"]:
        raise ValueError("Copied product inventory differs from recorded normalization")
    native.update(phylogeny=str(view), products={name: record(view / name)
                  for name in (*PRODUCTS, "species_tree.rooted.nwk")})
    pins = [manifest_pin, *proof["original"].values(), *proof["derived"].values(),
            *copied, manifest["original_completed_native_metadata"]]
    return native, projection, pins


def cell(repo, root, plan, report):
    require_completed(report)
    scheduler = terminal(report)
    original, parent, failed, parent_pins = parent_context(repo, plan["failure"])
    native, projection, native_pins = export_readback(root, report, original, parent)
    original_report = record(Path(original["output"]) / "results.json")
    control = {"report": failed["control"]["report"], "scores": failed["control"]["cell"]["scores"]}
    if (report.get("native_products") != native["products"]
            or report.get("native_summary") != native["summary"]
            or report.get("coverage") != {k: native[k] for k in ("genes", "species", "groups", "pairs")}
            or report.get("original_failed_producer") != original_report
            or report.get("control") != control
            or report.get("converter") != record(Path(plan["source_root"]) / "qfo_benchmark/og_to_pairwise.py")):
        raise ValueError("Recovered scores do not belong to authenticated native exports and matched control")
    scores, details, score_pins = readback_scores(repo, "QfO", report["scoring"], native, root,
                                               "ohmm_hillclimb_r1_" + plan["scoring_id"])
    vector = endpoint_vector("QfO", scores)
    if (scores != report["scores"] or report.get("deltas") != {
            key: value - control["scores"][key] for key, value in vector.items()}):
        raise ValueError("Recovered shared-QfO scores or matched-control deltas differ")
    time_pin = record(root / "recovery.time.txt")
    resources = verbose_time(Path(time_pin["path"]).read_text())
    checked = [report["execution_plan"], record(root / "results.json"), plan["failure"],
               original_report, time_pin, report["converter"], *parent_pins, *native_pins, *score_pins]
    for pin in checked:
        check(pin)
    compact = [report["execution_plan"], record(root / "results.json"), plan["failure"], original_report,
               time_pin, report["converter"], control["report"], *native_pins, *score_pins]
    return {"dataset": "QfO", "status": "independent_content_readback_verified", "scores": vector,
        "details": details, "products": native["products"], "coverage": report["coverage"],
        "checked": list({p["path"]: p for p in compact}.values()),
        "verification_record_count": len({p["path"] for p in checked}),
        "source_revision": REVISION, "matched_control": control,
        "normalization": report["normalization"], "independent_projection": projection,
        "resources": {"completed_native_stage_in_failed_export_job": failed["execution"],
            "derived_export_validation_and_scoring": {"scheduler": scheduler, "gnu_time": resources, "source": time_pin}},
        "execution_kind": "lossless_derived_export_scoring_of_original_canonical_shared_hit1_cap100",
        "scope": "Original failed producer and skipped attempt retained. No inference rerun or distinct relation/confidence change. Separate recovery costs are not uninterrupted full-pipeline timing. Duplicate cause remains unknown; FAS remains unseeded."}


def require_existing_control(methods, cell):
    controls = [m for m in methods if m["key"] == "canonical_collection:canonical"]
    scores = cell["matched_control"]["scores"]
    if (len(controls) != 1 or len(scores) != 6
            or {k: v for k, v in controls[0]["scores"].items() if k.startswith("QfO:")} != scores):
        raise ValueError("Existing canonical control differs from independently reverified matched control")

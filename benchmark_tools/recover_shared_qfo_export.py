"""Recover scoring through an explicit lossless view of a rejected native export."""

import argparse
import csv
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.continue_shared_qfo import REVISION, relocated_case
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.recover_cap400_scores import verify_helpers, verify_imported_converter
from benchmark_tools.recover_full_qfo_search import short_id
from benchmark_tools.report_additional_hillclimb import canonical_recovery_cell, check_plan, shared_native
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores
from benchmark_tools.run_full_search_hillclimb import runtime, verify_input_inventory
from benchmark_tools.run_membership_hillclimb import PRODUCTS, score_qfo, write_json
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.shared_qfo_export import normalize, verify_projection

PROTOCOL = "hillclimb_shared_qfo_export_recovery_20261016.json"
FAILURE = "benchmark_tools/results/shared_qfo_duplicate_export_failure_20261016.json"
SLURM = "recover_shared_qfo_export.slurm"
DUPLICATE = ("tr|A7T186|A7T186_NEMVE", "UP000001593_45351",
             "tr|Q7QC33|Q7QC33_ANOGA", "UP000007062_7165")
ROWS, DISTINCT = 5561115, 5561114
DUPLICATE_LINE = 3869801
SUCCESS = "shared_qfo_derived_export_all_six_scores_verified_pending_admission"


def require_failure(failure, report, scheduler):
    expected = {"JobIDRaw": "24247", "State": "FAILED", "ExitCode": "1:0",
                "AllocCPUS": "16", "NodeList": "bizon"}
    stream = failure.get("pair_and_confidence_stream_audit", {})
    if (scheduler != [expected] or failure.get("job_id") != "24247" or report.get("job_id") != "24247"
            or failure.get("status") != "shared_qfo_terminal_duplicate_export_failure_reproduced_with_original_provenance"
            or failure.get("scientific_revision") != REVISION
            or report.get("status") != "failed_or_unverified_shared_qfo_continuation"
            or report.get("error_type") != "ValueError"
            or report.get("error") != "Native pairs are not canonical, unique and sorted"
            or report.get("scores") != {} or report.get("execution", {}).get("returncode") != 0
            or stream.get("rows_each") != ROWS or stream.get("distinct_pairs") != DISTINCT
            or stream.get("adjacent_duplicate_rows_each") != 1
            or stream.get("duplicate_confidence_conflicts") != 0):
        raise ValueError("Require the exact terminal shared-QfO redundant-export failure")
    return expected


def verify_control(repo, original, report):
    pin = report["control"]["report"]
    check(pin)
    root = Path(original["recovery_root"])
    if Path(pin["path"]).resolve() != root / "results.json":
        raise ValueError("Matched canonical control belongs to another recovery")
    control_report = json.loads(Path(pin["path"]).read_text())
    plan, pins = check_plan(repo, root, control_report)
    if control_report.get("job_id") != "24246" or control_report["execution_plan"] not in original["inputs"]:
        raise ValueError("Matched canonical control differs from the original dependency")
    control = canonical_recovery_cell(repo, plan, control_report, root)
    for key in ("scores", "products", "coverage", "source_revision"):
        if control[key] != report["control"]["cell"][key]:
            raise ValueError("Independent matched-control readback differs")
    return [pin, *pins, *control["checked"]]


def parent_context(repo, failure_pin):
    check(failure_pin)
    failure = json.loads(Path(failure_pin["path"]).read_text())
    for pin in failure["artifacts"]:
        check(pin)
    root = repo / "benchmarks/results/shared_qfo_continuation_r1_v1_20261014"
    plan_pin, report_pin = record(root / "execution_plan.json"), record(root / "results.json")
    if plan_pin not in failure["artifacts"] or report_pin not in failure["artifacts"]:
        raise ValueError("Failure evidence belongs to another parent")
    original = json.loads(Path(plan_pin["path"]).read_text())
    report = json.loads(Path(report_pin["path"]).read_text())
    text = subprocess.check_output(["sacct", "-j", "24247", "-X", "-P", "-o",
        "JobIDRaw,State,ExitCode,AllocCPUS,NodeList"], text=True)
    require_failure(failure, report, list(csv.DictReader(io.StringIO(text), delimiter="|")))
    checked = [failure_pin, plan_pin, report_pin, *failure["artifacts"], *original["inputs"],
               *original["sources"], *original["runtime"]["files"]]
    for pin in checked:
        check(pin)
    parent = json.loads(Path(original["parent_plan"]["path"]).read_text())
    if (original["repo"] != str(repo) or original["output"] != str(root)
            or original["scientific_revision"] != REVISION
            or relocated_case(parent, root) != original["case"]
            or report["execution"]["command"] != original["case"]["command"]
            or runtime() != original["runtime"] or original["runtime"] != parent["runtime"]):
        raise ValueError("Original scientific identity, command or runtime changed")
    verify_input_inventory(original["case"])
    folder = Path(original["case"]["output"])
    phylo, output = folder / "output/orthohmm_phylogeny", folder / "output"
    metrics = json.loads((folder / "metrics.json").read_text())
    for item in metrics["harness"]["output_manifest"]:
        pin = record(output / item["path"])
        if (pin["bytes"], pin["sha256"]) != (item["bytes"], item["sha256"]):
            raise ValueError("Original completed output manifest differs")
        checked.append(pin)
    checked += [record(folder / "metrics.json"), record(phylo / "species_tree.rooted.nwk"),
        *[record(p) for p in sorted((output / "orthohmm_working_res/high_sensitivity_checkpoint").iterdir())]]
    checked += verify_control(repo, original, report)
    return original, parent, report, list({p["path"]: p for p in checked}.values())


def prepare(repo, source, output):
    if source == repo or output.exists():
        raise ValueError("Require separate frozen executor and fresh recovery namespace")
    failure_pin = record(repo / FAILURE)
    if record(source / FAILURE)["sha256"] != failure_pin["sha256"]:
        raise ValueError("Frozen failure evidence differs")
    original, _, _, _ = parent_context(repo, failure_pin)
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    pins = verify_helpers(source, original, revision, protocol=PROTOCOL, slurm=SLURM)
    identifier = short_id(repo, output, 1)
    if (repo / "qfo_benchmark/w" / ("hc1_" + identifier)).exists():
        raise FileExistsError("Scoring work namespace already exists")
    result = {"status": "prepared_shared_qfo_lossless_export_recovery", "repo": str(repo),
        "source_root": str(source), "source_revision": revision, "output": str(output),
        "failure": failure_pin, "sources": pins, "runtime": original["runtime"],
        "scientific_revision": REVISION, "scoring_id": identifier, "native_inference_rerun": False,
        "scope": "Explicit derived canonical-export scoring; original failed inference/export history retained, not uninterrupted successful pipeline timing."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def derive(original, parent, output):
    phylo = Path(original["case"]["output"]) / "output/orthohmm_phylogeny"
    view = output / "canonical_export"
    view.mkdir(exist_ok=False)
    proof = normalize(phylo / PRODUCTS[1], phylo / PRODUCTS[2], view / PRODUCTS[1], view / PRODUCTS[2],
                      DUPLICATE, ROWS, DISTINCT)
    if (proof["duplicate_confidence"] != "medium"
            or any(proof["skips"][name]["line"] != DUPLICATE_LINE for name in ("pairs", "confidence"))):
        raise ValueError("Original redundant-row location or confidence changed")
    projection = verify_projection(proof)
    native = shared_native(original["case"], parent, pair_export=proof)
    if ((native["genes"], native["species"], native["pairs"]) != (984137, 78, DISTINCT)
            or native["status"] != "derived_native_content_verified"):
        raise ValueError("Derived export lost complete original native coverage")
    copied = []
    for name in (PRODUCTS[0], "provenance_manifest.json", "reconciliation_summary.json", "species_tree.rooted.nwk"):
        before = record(phylo / name)
        shutil.copyfile(phylo / name, view / name)
        after = record(view / name)
        if (before["bytes"], before["sha256"]) != (after["bytes"], after["sha256"]):
            raise ValueError("Derived metadata or root groups differ from original")
        copied += [before, after]
    native.update(phylogeny=str(view), products={name: record(view / name)
                  for name in (*PRODUCTS, "species_tree.rooted.nwk")})
    manifest = {"status": "explicit_lossless_canonical_export_view_verified", "normalization": proof,
        "independent_projection": projection, "unchanged_copied_products": copied,
        "original_completed_native_metadata": record(Path(original["case"]["output"]) / "metrics.json"),
        "scope": "Copied native metadata describes the original inference, not a new execution. Only the redundant export row is omitted from each derived table.",
        "native_inference_rerun": False}
    write_json(output / "normalization.json", manifest)
    return native, manifest


def run(path):
    plan_pin = record(path)
    plan = json.loads(path.read_text())
    repo, source, output = (Path(plan[k]) for k in ("repo", "source_root", "output"))
    if (Path(__file__).resolve() != source / "benchmark_tools/recover_shared_qfo_export.py"
            or path.resolve() != output / "execution_plan.json" or (output / "results.json").exists()):
        raise ValueError("Wrong frozen executor or attempted recovery namespace")
    if (os.environ.get("SLURM_CPUS_PER_TASK"), os.environ.get("SLURM_MEM_PER_NODE"),
            os.environ.get("SLURMD_NODENAME"), os.environ.get("CUDA_VISIBLE_DEVICES")) != (
            "8", "65536", "bizon", ""):
        raise ValueError("Require CPU-only 8-CPU 64-GiB bizon allocation")
    result = {"status": "validating_failed_shared_qfo_for_lossless_export_recovery",
        "execution_plan": plan_pin, "job_id": os.environ.get("SLURM_JOB_ID"), "scores": {},
        "source_revision": plan["source_revision"], "scientific_revision": REVISION,
        "native_inference_rerun": False, "defaults_promoted": False, "goal_completed": False,
        "host_before": host_activity(), "scope": plan["scope"]}
    write_json(output / "results.json", result)
    checked = [plan_pin, plan["failure"], *plan["sources"], *plan["runtime"]["files"]]
    started = time.monotonic()
    try:
        for pin in checked:
            check(pin)
        original, parent, failed, parent_pins = parent_context(repo, plan["failure"])
        checked += parent_pins
        if (plan.get("native_inference_rerun") is not False or plan["scientific_revision"] != REVISION
                or plan["runtime"] != original["runtime"] or short_id(repo, output, 1) != plan["scoring_id"]
                or verify_helpers(source, original, plan["source_revision"], protocol=PROTOCOL, slurm=SLURM) != plan["sources"]):
            raise ValueError("Recovery changed original identity, executor or namespace")
        result["converter"] = verify_imported_converter(source)
        result["status"] = "deriving_and_validating_lossless_native_export"
        write_json(output / "results.json", result)
        native, manifest = derive(original, parent, output)
        checked += [record(output / "normalization.json"), *manifest["unchanged_copied_products"],
                    *manifest["normalization"]["derived"].values()]
        result.update(status="scoring_derived_shared_qfo_relation_set", normalization=record(output / "normalization.json"),
            coverage={k: native[k] for k in ("genes", "species", "groups", "pairs")},
            native_products=native["products"], native_summary=native["summary"],
            original_failed_producer=record(Path(original["output"]) / "results.json"),
            control={"report": failed["control"]["report"], "scores": failed["control"]["cell"]["scores"]},
            export_verification_wall_s=time.monotonic() - started)
        write_json(output / "results.json", result)
        scoring_started = time.monotonic()
        scored = score_qfo(repo, Path(native["phylogeny"]), output, plan["scoring_id"])
        scores, details, pins = readback_scores(repo, "QfO", scored, native, output,
                                               "ohmm_hillclimb_r1_" + plan["scoring_id"])
        vector = endpoint_vector("QfO", scores)
        checked += pins
        result.update(status=SUCCESS, scores=scores, details=details, scoring=scored,
            deltas={k: value - result["control"]["scores"][k] for k, value in vector.items()},
            scoring_and_readback_wall_s=time.monotonic() - scoring_started)
    except Exception as error:
        result.update(status="failed_shared_qfo_lossless_export_recovery", error_type=type(error).__name__, error=str(error))
    finally:
        try:
            for pin in checked:
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Original runtime changed during export recovery")
            result["postflight_pins_verified"] = True
        except Exception as error:
            result.update(status="failed_shared_qfo_lossless_export_recovery", postflight_error=str(error))
        result.update(total_recovery_wall_s=time.monotonic() - started, host_after=host_activity())
        write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    preparation = sub.add_parser("prepare")
    for name in ("repo", "source", "output"):
        preparation.add_argument("--" + name, required=True, type=Path)
    execution = sub.add_parser("run")
    execution.add_argument("--execution-plan", required=True, type=Path)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve())
              if args.operation == "prepare" else run(args.execution_plan.resolve()))
    print(result["status"])
    raise SystemExit(0 if result["status"] in ("prepared_shared_qfo_lossless_export_recovery", SUCCESS) else 1)

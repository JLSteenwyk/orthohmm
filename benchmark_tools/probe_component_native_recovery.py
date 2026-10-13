"""Validate completed component inference in a clean process after parent failure."""

import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json


def validate_native_receipt(metrics, command, revision, time_text):
    if metrics.get("status") != "complete" or metrics.get("command") != command:
        raise ValueError("Retained native execution did not complete the planned command")
    if metrics.get("git") != {"commit": revision, "dirty": False}:
        raise ValueError("Native execution source differs from frozen source")
    if [line.strip() for line in time_text.splitlines() if "Exit status:" in line] != ["Exit status: 0"]:
        raise ValueError("Native GNU time has no unambiguous zero exit status")


def run(path, output):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Require new recovery namespace and CPU-only runtime")
    plan_pin = record(path)
    plan = json.loads(path.read_text())
    source, folder = Path(plan["source_root"]), Path(plan["output"])
    producer_path = folder / "results.json"
    producer_pin = record(producer_path)
    producer = json.loads(producer_path.read_text())
    if plan["dataset"] != "QfO" or producer["execution_plan"] != plan_pin or producer["source_revision"] != plan["source_revision"]:
        raise ValueError("Mismatched original component execution plan")
    accounting = subprocess.check_output(["sacct", "-j", producer["job_id"], "-X",
        "--format=JobID,State,ExitCode", "--parsable2"], text=True)
    jobs = list(csv.DictReader(accounting.splitlines(), delimiter="|"))
    if len(jobs) != 1 or jobs[0]["JobID"] != producer["job_id"] or jobs[0]["State"] != "FAILED" or jobs[0]["ExitCode"] != "0:11":
        raise ValueError("Original parent is not authoritatively terminal with signal 11")
    case = plan["case"]
    native_folder = Path(case["output"])
    metrics_path, time_path = native_folder / "metrics.json", native_folder / "native.time.txt"
    metrics_pin, time_pin = record(metrics_path), record(time_path)
    metrics = json.loads(metrics_path.read_text())
    validate_native_receipt(metrics, case["command"], plan["source_revision"], time_path.read_text())
    checked = [plan_pin, producer_pin, metrics_pin, time_pin, *plan["sources"], *plan["runtime"]["files"]]
    checked.extend(producer["structural_controls"]["components500"]["candidate"][key]
                   for key in ("candidate_partition", "membership_constraints"))
    for pin in checked:
        check(pin)
    output.mkdir(parents=True)
    report = {"status": "clean_native_readback_running", "scheduler": jobs[0], "checked": checked,
              "probe": record(Path(__file__)), "scores": {}, "inference_rerun": False,
              "defaults_promoted": False, "goal_completed": False}
    write_json(output / "report.json", report)
    try:
        sys.path.insert(0, str(source))
        spec = importlib.util.spec_from_file_location("component_original_native_reader", source / "benchmark_tools/report_hillclimb.py")
        reader = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reader)
        from benchmark_tools.run_membership_hillclimb import validate_partition
        from benchmark_tools.run_copy_branch_ablation import require_markers
        phylo = native_folder / "output/orthohmm_phylogeny"
        candidate = producer["structural_controls"]["components500"]["candidate"]["candidate_partition"]
        coverage = validate_partition(Path(candidate["path"]), phylo.parent / "orthohmm_working_res/orthohmm_edges_clustered.txt")
        row = {"outputs": [record(phylo / name) for name in reader.PRODUCTS], "coverage": coverage}
        native = reader.replay_native(row, case)
        require_markers(native["summary"])
        expected_summary = {**metrics["counts"], "schema_version": 1}
        if native["genes"] != 984137 or native["species"] != 78 or native["summary"] != expected_summary:
            raise ValueError("Native receipt and reconstructed full-dataset counts disagree")
        for pin in [*checked, *native["checked"]]:
            check(pin)
        report.update(status="clean_original_native_readback_verified_scoring_still_pending", native=native,
                      coverage=coverage, source_revision=plan["source_revision"])
    except Exception as error:
        report.update(status="failed_clean_native_readback", error_type=type(error).__name__, error=str(error))
    write_json(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution-plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.execution_plan.resolve(), args.output.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"].startswith("clean_original_") else 1)

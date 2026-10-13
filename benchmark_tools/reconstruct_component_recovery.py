"""Independently reconstruct frozen component candidates in an isolated process."""

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys


def compare_reconstruction(previous, rows, diagnostics, check):
    if set(rows) != {"control", "components500"}:
        raise ValueError("Incomplete reconstructed component arm inventory")
    canonical = json.loads(json.dumps(diagnostics["ownership"]))
    if (canonical != previous["ownership"] or diagnostics["component_counts"] != previous["component_counts"]):
        raise ValueError("Reconstructed ownership or component activation changed")
    checked = []
    for label in ("control", "components500"):
        old, new = previous["structural_controls"][label], rows[label]
        for key in ("seed", "copy_split_trace"):
            check(old[key])
            check(new[key])
            if old[key]["sha256"] != new[key]["sha256"]:
                raise ValueError("Reconstructed seed or executed copy-split trace changed")
            checked.extend((old[key], new[key]))
        for key in ("seed_partition", "candidate_partition", "membership_constraints"):
            a, b = old["candidate"][key], new["candidate"][key]
            check(a)
            check(b)
            if a["sha256"] != b["sha256"]:
                raise ValueError("Reconstructed native candidate or constraints changed")
            checked.extend((a, b))
    return checked


def run(path, output):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Require a new CPU-only reconstruction namespace")
    plan = json.loads(path.read_text())
    source = Path(plan["source_root"])
    sys.path.insert(0, str(source))
    from benchmark_tools import run_large_copy_components as runner
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.run_membership_hillclimb import write_json
    if Path(runner.__file__).resolve() != source / "benchmark_tools/run_large_copy_components.py":
        raise ValueError("Reconstruction loaded the wrong original scientific source")
    producer_path = Path(plan["output"]) / "results.json"
    previous = json.loads(producer_path.read_text())
    if previous["execution_plan"] != record(path) or previous["source_revision"] != plan["source_revision"]:
        raise ValueError("Reconstruction belongs to another original execution")
    accounting = subprocess.check_output(["sacct", "-j", previous["job_id"], "-X",
        "--format=JobID,State,ExitCode", "--parsable2"], text=True)
    jobs = list(csv.DictReader(accounting.splitlines(), delimiter="|"))
    if (len(jobs) != 1 or jobs[0]["JobID"] != previous["job_id"]
            or jobs[0]["State"] != "FAILED" or jobs[0]["ExitCode"] != "0:11"):
        raise ValueError("Original component parent is not terminal with signal 11")
    checked = [record(path), record(producer_path), record(Path(__file__)),
               record(Path(__file__).with_name("hillclimb_component_candidate_recovery_20261013.json"))]
    runner.verify(plan, checked)
    output.mkdir(parents=True)
    report = {"status": "isolated_component_reconstruction_running", "job_id": os.environ.get("SLURM_JOB_ID"),
        "original_execution_plan": record(path), "original_producer": record(producer_path), "scheduler": jobs[0],
        "source_revision": plan["source_revision"], "executor": record(Path(__file__)), "scores": {},
        "native_inference_rerun": False, "defaults_promoted": False, "goal_completed": False}
    write_json(output / "report.json", report)
    try:
        rows, diagnostics = runner.construct(plan, output)
        checked.extend(compare_reconstruction(previous, rows, diagnostics, check))
        if diagnostics["component_counts"] != plan["expected_component_counts"]:
            raise ValueError("Component activation differs from prespecified feasibility counts")
        runner.verify(plan, checked)
        report.update(status="isolated_component_candidates_exactly_reconstructed", structural_controls=rows,
                      component_counts=diagnostics["component_counts"], ownership=json.loads(json.dumps(diagnostics["ownership"])),
                      checked=checked, inputs_reverified_unchanged=True)
    except Exception as error:
        report.update(status="failed_isolated_component_reconstruction", error_type=type(error).__name__, error=str(error))
    write_json(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.original_plan.resolve(), args.output.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"] == "isolated_component_candidates_exactly_reconstructed" else 1)

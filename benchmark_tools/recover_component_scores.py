"""Recover QfO scoring from independently verified, retained component outputs."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def require_completed_probe(text, job_id):
    rows = list(csv.DictReader(text.splitlines(), delimiter="|"))
    if (len(rows) != 1 or rows[0]["JobID"] != str(job_id)
            or rows[0]["State"] != "COMPLETED" or rows[0]["ExitCode"] != "0:0"
            or rows[0]["AllocCPUS"] != "1"):
        raise ValueError("Require successful terminal one-CPU native verification probe")
    return rows[0]


def original_helpers(source):
    sys.path.insert(0, str(source))
    from benchmark_tools import run_large_copy_components as runner
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.run_membership_hillclimb import score_qfo, write_json
    if Path(runner.__file__).resolve() != source / "benchmark_tools/run_large_copy_components.py":
        raise ValueError("Recovery loaded the wrong frozen scientific helpers")
    return runner, check, record, score_qfo, write_json


def prepare(original_path, probe_path, output):
    original = json.loads(original_path.read_text())
    runner, check, record, _, write_json = original_helpers(Path(original["source_root"]))
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Require a fresh recovery namespace and CPU-only runtime")
    probe = json.loads(probe_path.read_text())
    if probe["status"] != "clean_original_native_readback_verified_scoring_still_pending" or probe["inference_rerun"]:
        raise ValueError("Native products have not passed clean original readback")
    accounting = subprocess.check_output(["sacct", "-j", "24219", "-X",
        "--format=JobID,State,ExitCode,AllocCPUS", "--parsable2"], text=True)
    scheduler = require_completed_probe(accounting, "24219")
    pin = record(original_path)
    if pin not in probe["checked"] or probe["source_revision"] != original["source_revision"]:
        raise ValueError("Clean native verification belongs to another execution plan")
    if probe["native"]["genes"] != 984137 or probe["native"]["species"] != 78:
        raise ValueError("Recovery lost full QfO gene or species coverage")
    runner.verify(original, [*probe["checked"], *probe["native"]["checked"], record(probe_path)])
    index = "cpr" + hashlib.sha256(str(output).encode()).hexdigest()[:10]
    from benchmark_tools.run_qfo_recovered_assessment import require_darwin_work_path
    repo = Path(original["repo"])
    work = repo / "qfo_benchmark/w" / ("hc1_" + index)
    require_darwin_work_path(work)
    if work.exists():
        raise FileExistsError("QfO recovery work namespace already exists")
    executor = Path(__file__).resolve()
    executor_root = executor.parent.parent
    revision = subprocess.check_output(["git", "-C", str(executor_root), "rev-parse", "HEAD"], text=True).strip()
    producer_path = Path(original["output"]) / "results.json"
    plan = {"status": "prepared_scoring_only_after_clean_native_verification", "original_execution_plan": pin,
        "original_failed_producer": record(producer_path), "original_source_root": original["source_root"],
        "original_source_revision": original["source_revision"], "native_probe": record(probe_path),
        "native_probe_scheduler": scheduler, "executor": record(executor), "executor_root": str(executor_root),
        "executor_revision": revision, "protocol": record(executor.with_name("hillclimb_component_scoring_recovery_20261013.json")),
        "output": str(output), "scoring_id": index, "work": str(work), "dataset": "QfO",
        "native_inference_rerun": False, "defaults_promoted": False, "goal_completed": False}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", plan)
    return plan


def run(path):
    plan = json.loads(path.read_text())
    original = json.loads(Path(plan["original_execution_plan"]["path"]).read_text())
    runner, check, record, score_qfo, write_json = original_helpers(Path(plan["original_source_root"]))
    output = Path(plan["output"])
    if (output / "results.json").exists() or record(Path(__file__)) != plan["executor"]:
        raise ValueError("Require untouched recovery execution and pinned executor")
    probe = json.loads(Path(plan["native_probe"]["path"]).read_text())
    native = probe["native"]
    checked = [record(path), plan["original_execution_plan"], plan["original_failed_producer"],
               plan["native_probe"], plan["executor"], plan["protocol"], *probe["checked"], *native["checked"]]
    for pin in checked:
        check(pin)
    runner.verify(original, checked)
    report = {"status": "scoring_recovery_running", "execution_plan": record(path),
              "source_revision": plan["original_source_revision"], "job_id": os.environ.get("SLURM_JOB_ID"),
              "scores": {}, "native_inference_rerun": False, "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", report)
    started = time.monotonic()
    try:
        scored = score_qfo(Path(original["repo"]), Path(native["phylogeny"]), output, plan["scoring_id"])
        from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores
        scores, details, files = readback_scores(Path(original["repo"]), "QfO", scored, native, output,
            "ohmm_hillclimb_r1_" + plan["scoring_id"])
        runner.verify(original, [*checked, *files])
        report.update(status="scored_recovery_pending_independent_admission", scores=endpoint_vector("QfO", scores),
            details=details, scoring=scored, native_products=native["products"], native_summary=native["summary"],
            coverage={k: native[k] for k in ("genes", "species", "groups", "pairs")},
            scope="QfO scoring-only recovery on retained candidate-specific own-tree inference; no native rerun or fresh pipeline timing.")
    except Exception as error:
        report.update(status="failed_scoring_only_recovery", error_type=type(error).__name__, error=str(error))
    report["scoring_and_verification_wall_s"] = time.monotonic() - started
    write_json(output / "results.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    preparation = sub.add_parser("prepare")
    preparation.add_argument("--original-plan", type=Path, required=True)
    preparation.add_argument("--native-probe", type=Path, required=True)
    preparation.add_argument("--output", type=Path, required=True)
    execution = sub.add_parser("run")
    execution.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    if args.operation == "prepare":
        result = prepare(args.original_plan.resolve(), args.native_probe.resolve(), args.output.resolve())
    else:
        result = run(args.execution_plan.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"] in ("prepared_scoring_only_after_clean_native_verification",
                                              "scored_recovery_pending_independent_admission") else 1)

"""Finish the original unrun cap400 arm after authenticated native recovery."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.recover_full_qfo_search import (
    parent_cases, planned_case, verified_native, independent_scores,
)
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
from benchmark_tools.run_full_search_hillclimb import checked_records, runtime, verify_input_inventory
from benchmark_tools.run_membership_hillclimb import score_qfo, write_json
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import sources, timed_stage

PROTOCOL = "hillclimb_qfo_cap400_continuation_20261014.json"


def require_recovered_control(report):
    if report["status"] != "cached_native_control_and_all_six_scores_verified":
        raise ValueError("Recovered original control and all six scores must be independently verified")
    endpoint_vector("QfO", report["scores"])


def prepare(repo, source, output):
    if output.exists() or source == repo:
        raise ValueError("Require fresh continuation namespace and frozen runner")
    parent_path = repo / "benchmarks/results/full_search_hillclimb_r1_v3_20261012/execution_plan.json"
    parent = json.loads(parent_path.read_text())
    original = parent_cases(parent)[1]
    continuation = planned_case(repo, output, original, 1)
    control_root = repo / "benchmarks/results/qfo_native_recovery_r1_v1_20261014"
    inputs = [record(source / "benchmark_tools" / PROTOCOL), record(parent_path),
              record(control_root / "execution_plan.json"), *checked_records(parent, original)]
    source_pins = sources(source) + [record(source / "benchmark_tools/continue_full_qfo_cap400.slurm")]
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    for pin in source_pins:
        relative = Path(pin["path"]).relative_to(source)
        if relative.suffix == ".so": continue
        blob = subprocess.check_output(["git", "show", revision + ":" + str(relative)], cwd=repo)
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Continuation runner differs from committed source")
    for pin in (*inputs, *source_pins): check(pin)
    if runtime() != parent["runtime"]:
        raise ValueError("Original scientific runtime changed")
    result = {"status": "prepared_waiting_for_recovered_control", "repo": str(repo), "source_root": str(source),
              "source_revision": revision, "scientific_source": parent["source_root"],
              "scientific_revision": parent["source_revision"], "output": str(output), "case": continuation,
              "control_root": str(control_root), "inputs": inputs, "sources": source_pins,
              "runtime": parent["runtime"], "tool_path": parent["tool_path"], "dependency_job": "24244",
              "scope": "Original fresh cap400 continuation; failed full control and cached recovery costs remain separate."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def run(path):
    plan = json.loads(path.read_text())
    source, science, repo, output = (Path(plan[key]) for key in ("source_root", "scientific_source", "repo", "output"))
    if Path(__file__).resolve() != source / "benchmark_tools/continue_full_qfo_cap400.py" or (output / "results.json").exists():
        raise ValueError("Require frozen runner and preserve previous attempts")
    result = {"status": "validating_control", "execution_plan": record(path), "job_id": os.environ.get("SLURM_JOB_ID"),
              "scope": plan["scope"], "scores": {}, "defaults_promoted": False, "goal_completed": False,
              "host_before": host_activity()}
    write_json(output / "results.json", result)
    try:
        accounting = subprocess.check_output(["sacct", "-j", plan["dependency_job"], "-X", "-n", "-P",
                                             "-o", "JobID,State%30,ExitCode"], text=True)
        if accounting.strip() != plan["dependency_job"] + "|COMPLETED|0:0":
            raise ValueError("Recovery scheduler is not terminal successful")
        control_root = Path(plan["control_root"])
        control_path = control_root / "results.json"
        control = json.loads(control_path.read_text())
        require_recovered_control(control)
        check(control["execution_plan"])
        control_plan = json.loads(Path(control["execution_plan"]["path"]).read_text())
        native_control = replay_native(control, {"command": control_plan["command"]})
        scoring = json.loads((control_root / "scoring_preflight.json").read_text())
        control_scores, control_details, pins = readback_scores(repo, "QfO", control, native_control, control_root,
                                                               scoring["stage"]["participant"])
        result["control"] = {"report": record(control_path), "scores": control_scores, "details": control_details,
                              "checked": [*native_control["checked"], *pins]}
        checked = [record(path), record(control_path), *plan["inputs"], *plan["sources"],
                   *plan["runtime"]["files"], *result["control"]["checked"]]
        for pin in checked: check(pin)
        if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
            raise ValueError("Changed runtime or CUDA visibility")
        case = plan["case"]["fresh_case"]
        verify_input_inventory(case)
        native_folder = Path(case["output"])
        native_folder.mkdir(parents=True, exist_ok=False)
        result["status"] = "fresh_original_cap400_running"
        write_json(output / "results.json", result)
        env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(science), "PYTHONHASHSEED": "0",
               "CUDA_VISIBLE_DEVICES": "", "PYTHONFAULTHANDLER": "1", "OMP_NUM_THREADS": "1",
               "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "32"}
        result["execution"] = timed_stage(case["command"], native_folder, "native", science, env)
        if result["execution"]["returncode"]:
            raise RuntimeError("Original fresh cap400 failed; retain artifacts for separate recovery")
        native = verified_native(case)
        result["native"] = native
        result["status"] = "native_verified_scoring_pending"
        write_json(output / "results.json", result)
        scoring_folder = Path(plan["case"]["scoring_output"])
        scoring_folder.mkdir(parents=True, exist_ok=False)
        result.update(score_qfo(repo, Path(native["phylogeny"]), scoring_folder, plan["case"]["scoring_id"]))
        result["independent_scores"] = independent_scores(repo, result, native, scoring_folder, plan["case"]["scoring_id"])
        result["deltas"] = {name: value - control_scores[name] for name, value in result["scores"].items()}
        for pin in (*checked, *native["checked"], *result["independent_scores"]["checked"]): check(pin)
        if runtime() != plan["runtime"]:
            raise ValueError("Runtime changed during continuation")
        result["status"] = "original_cap400_native_and_all_six_scores_verified"
    except Exception as error:
        result.update(status="failed_or_unverified_control_continuation", error_type=type(error).__name__, error=str(error))
    result["host_after"] = host_activity()
    write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    preparation = sub.add_parser("prepare")
    for name in ("repo", "source", "output"): preparation.add_argument("--" + name, type=Path, required=True)
    worker = sub.add_parser("run")
    worker.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve()) if args.command == "prepare" else run(args.execution_plan.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"] in {"prepared_waiting_for_recovered_control", "original_cap400_native_and_all_six_scores_verified"} else 1)

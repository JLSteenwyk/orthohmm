"""Recover scoring of completed fresh QfO searches without repeating valid inference."""

from __future__ import annotations

import argparse
from copy import deepcopy
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_full_search_hillclimb import (
    checked_records, native_command, runtime, validate_native, verify_input_inventory,
)
from benchmark_tools.run_held_tree_hillclimb import source_records
from benchmark_tools.run_membership_hillclimb import argument, score_qfo, write_json
from benchmark_tools.run_qfo_recovered_assessment import require_darwin_work_path
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.validate_profile_runtime import require_profile_runtime

REVISION = "2ccc411520d9e63678178a50207ba1fca87bfbb3"
LABELS = ("control100", "cap400")
PROTOCOL = "hillclimb_full_qfo_scoring_recovery_20261013.json"
PATH_ERROR = "Darwin work path exceeds validated length limit"


def parent_cases(parent):
    cases = [c for c in parent["cases"] if c["dataset"] == "QfO"]
    if (parent.get("source_revision") != REVISION or type(parent.get("cpu")) is not int
            or parent["cpu"] != 32 or tuple(c["label"] for c in cases) != LABELS
            or [c["search_max_candidates_per_query"] for c in cases] != [100, 400]):
        raise ValueError("Require the frozen original 32-CPU QfO cap trial")
    source = Path(parent["source_root"])
    for case in cases:
        if len(case["inputs"]) != 78:
            raise ValueError("Require the complete 78-proteome QfO inventory")
        tools = [argument(case["command"], flag) for flag in ("--aligner", "--tree-builder")]
        expected = native_command(source, Path(case["fasta"]), Path(case["output"]), tools,
                                  case["search_max_candidates_per_query"], 32)
        if case["command"] != expected:
            raise ValueError("Parent native command differs from the prespecified trial")
    if cases[0]["inputs"] != cases[1]["inputs"] or cases[0]["fasta"] != cases[1]["fasta"]:
        raise ValueError("Parent arms use different input inventories")
    return cases


def terminal_scheduler(accounting, job):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    matching = [r for r in rows if r.get("JobIDRaw") == str(job)]
    if len(matching) != 1:
        raise ValueError("Parent job is not uniquely present in scheduler accounting")
    row = matching[0]
    if ((row.get("State"), row.get("ExitCode")) not in {("FAILED", "1:0"), ("COMPLETED", "0:0")}
            or row.get("NodeList") != "bizon" or row.get("AllocCPUS") != "32"):
        raise ValueError("Require terminal parent with expected exit, host and 32 CPUs")
    return row


def recovery_actions(report, parent_identity, job, scheduler):
    if (report.get("job_id") != str(job) or report.get("execution_plan") != parent_identity
            or report.get("source_revision") != REVISION
            or report.get("status") not in {"finished_all_arms_scored", "finished_with_explicit_gaps"}
            or tuple(r.get("label") for r in report.get("arms", [])) != LABELS):
        raise ValueError("Require the exact terminal parent report and ordered arm inventory")
    actions = []
    for index, row in enumerate(report["arms"]):
        if row.get("inputs_sources_runtime_reverified") is not True:
            raise ValueError("Parent postflight input/source/runtime gate failed")
        status = row.get("status")
        if status == "not_run_failed_control" and index == 1 and actions == ["score_existing"]:
            if row.get("scores") or "native" in row or "exit_code" in row:
                raise ValueError("Skipped variant claims inference or scores")
            actions.append("run_missing")
            continue
        if (row.get("exit_code") != 0 or row.get("cpu_backend_probe", {}).get("passed") is not True
                or row.get("native", {}).get("status") != "fresh_native_content_verified"):
            raise ValueError("Parent arm has no verified successful fresh native inference")
        if status == "scored_development_result":
            actions.append("reuse_scored")
        elif (status == "failed" and row.get("error_type") == "ValueError"
              and row.get("error") == PATH_ERROR and not row.get("scores")):
            actions.append("score_existing")
        else:
            raise ValueError("Parent failure is not the specific recoverable scoring-path error")
    expected = "COMPLETED" if actions == ["reuse_scored", "reuse_scored"] else "FAILED"
    if scheduler["State"] != expected or report["status"] != (
            "finished_all_arms_scored" if expected == "COMPLETED" else "finished_with_explicit_gaps"):
        raise ValueError("Scheduler and terminal report disagree")
    return actions


def short_id(repo, output, index):
    if type(index) is not int or index not in (0, 1):
        raise ValueError("Require a planned arm index")
    name = "fq_" + hashlib.sha256(str(output).encode()).hexdigest()[:16] + f"_{index}"
    require_darwin_work_path(repo / "qfo_benchmark/w" / ("hc1_" + name))
    return name


def planned_case(repo, output, original, index):
    folder = output / original["label"]
    native = deepcopy(original)
    native["output"] = str(folder / "native")
    command = list(original["command"])
    command[3:5] = [str(folder / "native/output"), str(folder / "native/metrics.json")]
    native["command"] = command
    return {"label": original["label"], "original_case": original, "fresh_case": native,
            "scoring_output": str(folder / "scoring"), "scoring_id": short_id(repo, output, index)}


def prepare(repo, source, output, parent_path, parent_job):
    if output.exists() or not str(parent_job).isdigit() or int(parent_job) <= 0:
        raise ValueError("Require a fresh namespace and explicit positive parent job")
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools"], check=True)
    parent = json.loads(parent_path.read_text())
    cases = parent_cases(parent)
    if Path(parent["repo"]).resolve() != repo:
        raise ValueError("Parent belongs to another repository")
    for item in checked_records(parent):
        check(item)
    if runtime() != parent["runtime"]:
        raise ValueError("Parent runtime changed")
    require_profile_runtime(Path(parent["source_root"]))
    planned = [planned_case(repo, output, original, index) for index, original in enumerate(cases)]
    result = {"status": "prepared_waiting_for_terminal_parent", "repo": str(repo), "source_root": str(source),
        "output": str(output), "cpu": 32, "parent_plan": record(parent_path), "parent_job": str(parent_job),
        "parent_results": str(Path(parent["output"]) / "QfO/results.json"),
        "protocol": record(source / "benchmark_tools" / PROTOCOL), "cases": planned,
        "sources": [*source_records(source), record(source / "benchmark_tools/recover_full_qfo_search.slurm")],
        "source_revision": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "runtime": parent["runtime"], "tool_path": parent["tool_path"],
        "scope": "Original frozen fresh native search; completed outputs scored without repeating inference. Shared-host timings remain separated."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def verify_inputs(plan, parent, identities=()):
    for item in [plan["parent_plan"], plan["protocol"], *plan["sources"], *identities,
                 *checked_records(parent)]:
        check(item)
    for case in parent_cases(parent):
        verify_input_inventory(case)
    if runtime() != plan["runtime"] or plan["runtime"] != parent["runtime"]:
        raise ValueError("Frozen scientific runtime changed")
    source = Path(parent["source_root"])
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if revision != REVISION:
        raise ValueError("Original native source revision changed")
    require_profile_runtime(source)


def verified_native(case):
    native = validate_native(case, 32)
    if native["genes"] != 984137 or native["species"] != 78:
        raise ValueError("Wrong corrected QfO universe")
    harness = native["metrics"]["harness"]
    if harness.get("git_commit") != REVISION or harness.get("exit_code") != 0:
        raise ValueError("Wrong completed native source revision")
    source = Path(case["command"][1]).parent.parent
    expected_sources = [record(p) for p in sorted((source / "orthohmm").rglob("*.py"))]
    expected_sources.append(record(source / "benchmark_tools/benchmark_production.py"))
    actual = {r["path"]: (r["bytes"], r["sha256"]) for r in harness["source_manifest"]}
    expected = {str(Path(r["path"]).relative_to(source)): (r["bytes"], r["sha256"]) for r in expected_sources}
    inputs = {Path(r["path"]).name: (r["bytes"], r["sha256"]) for r in case["inputs"]}
    if (len(actual) != len(harness["source_manifest"]) or actual != expected
            or len(harness["input_manifest"]) != 78
            or {r["path"]: (r["bytes"], r["sha256"]) for r in harness["input_manifest"]} != inputs):
        raise ValueError("Native harness source/input manifests differ from frozen trial")
    folder = Path(case["output"])
    checkpoint = folder / "output/orthohmm_working_res/high_sensitivity_checkpoint"
    checked = [*native["products"].values(), record(folder / "metrics.json"),
        record(folder / "output/orthohmm_phylogeny/provenance_manifest.json"),
        record(folder / "output/orthohmm_phylogeny/reconciliation_summary.json"),
        record(folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"),
        record(folder / "output/orthohmm_working_res/phylogeny_candidate_superfamilies.txt"),
        *(record(p) for p in sorted(checkpoint.iterdir()))]
    return {**native, "phylogeny": str(folder / "output/orthohmm_phylogeny"),
            "prediction": str(folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"), "checked": checked}


def independent_scores(repo, row, native, folder, scoring_id):
    from benchmark_tools.report_hillclimb import readback_scores
    scores, details, checked = readback_scores(repo, "QfO", row, native, folder,
                                             "ohmm_hillclimb_r1_" + scoring_id)
    for item in [*native["checked"], *checked]:
        check(item)
    return {"status": "independent_content_readback_verified", "scores": scores,
            "details": details, "checked": checked}


def resolve(plan, parent, output):
    if tuple(c.get("label") for c in plan["cases"]) != LABELS:
        raise ValueError("Require both ordered recovery cases")
    accounting = subprocess.check_output(["sacct", "-j", plan["parent_job"], "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = terminal_scheduler(accounting, plan["parent_job"])
    path = Path(plan["parent_results"])
    identity = record(path)
    report = json.loads(path.read_text())
    actions = recovery_actions(report, plan["parent_plan"], plan["parent_job"], scheduler)
    originals = parent_cases(parent)
    resolved, checked = [], [identity]
    for index, (row, case, planned, action) in enumerate(zip(report["arms"], originals, plan["cases"], actions)):
        if planned != planned_case(Path(plan["repo"]), output, case, index) or planned["label"] != row["label"]:
            raise ValueError("Recovery changed the original case inventory")
        if action == "run_missing":
            if Path(case["output"]).exists():
                raise ValueError("Skipped variant unexpectedly has native output; do not repeat it")
            native = None
        else:
            native = verified_native(case)
            if native["products"] != row["native"]["products"]:
                raise ValueError("Completed native products differ from parent receipt")
            checked.extend(native["checked"])
            if action == "score_existing":
                old_work = Path(plan["repo"]) / "qfo_benchmark/w" / ("hc1_" + case["scoring_id"])
                try:
                    require_darwin_work_path(old_work)
                except ValueError as error:
                    if str(error) != PATH_ERROR:
                        raise
                else:
                    raise ValueError("Reported path-length failure does not reproduce")
            else:
                readback = independent_scores(Path(plan["repo"]), row, native, Path(case["output"]), case["scoring_id"])
                checked.extend(readback["checked"])
        resolved.append({"label": case["label"], "action": action, "native": native, "parent_row": row})
    snapshot = output / "parent_results_snapshot.json"
    write_json(snapshot, report)
    checked.append(record(snapshot))
    for item in checked:
        check(item)
    receipt = {"status": "terminal_parent_native_outputs_resolved", "parent": identity,
               "snapshot": record(snapshot), "scheduler": scheduler, "checked": checked,
               "actions": [{"label": r["label"], "action": r["action"]} for r in resolved]}
    path = output / "resolved_parent_inputs.json"
    write_json(path, receipt)
    return resolved, [record(path), *checked]


def run_missing(case, parent, env):
    original = next(c for c in parent_cases(parent) if c["label"] == "cap400")
    expected = list(original["command"])
    folder = Path(case["output"])
    expected[3:5] = [str(folder / "output"), str(folder / "metrics.json")]
    if case != {**original, "output": str(folder), "command": expected}:
        raise ValueError("Only missing cap400 may run, with unchanged frozen native settings")
    folder.mkdir(parents=True, exist_ok=False)
    source = Path(parent["source_root"])
    subprocess.check_output([sys.executable, "-c",
        "from orthohmm.search import engine; from orthohmm.search.prefilter import is_c_prefilter_available; "
        "assert not engine.is_cuda_available() and engine.is_c_available() and is_c_prefilter_available()"],
        cwd=source, env=env, text=True)
    with (folder / "inference.log").open("x") as log:
        subprocess.run(case["command"], cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    return verified_native(case)


def run(path):
    plan = json.loads(path.read_text())
    source, output, repo = Path(plan["source_root"]), Path(plan["output"]), Path(plan["repo"])
    if Path(__file__).resolve() != source / "benchmark_tools/recover_full_qfo_search.py":
        raise ValueError("Execute the frozen recovery worker")
    if (output / "results.json").exists():
        raise FileExistsError("Recovery already started; inspect its live/terminal handle instead")
    report = {"status": "running", "job_id": os.environ.get("SLURM_JOB_ID"), "execution_plan": record(path),
        "source_revision": plan["source_revision"], "arms": [], "scope": plan["scope"], "defaults_promoted": False}
    write_json(output / "results.json", report)
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": "", "CUDA_VISIBLE_DEVICES": "",
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        "NUMBA_NUM_THREADS": "32", "PYTHONHASHSEED": "0"}
    identities = [report["execution_plan"]]
    try:
        parent = json.loads(Path(plan["parent_plan"]["path"]).read_text())
        env["PYTHONPATH"] = parent["source_root"]
        verify_inputs(plan, parent, identities)
        resolved, pins = resolve(plan, parent, output)
        identities.extend(pins)
        report["resolved_parent_inputs"] = pins[0]
        for old, case in zip(resolved, plan["cases"]):
            row = {"label": old["label"], "action": old["action"], "status": "running", "scores": {}}
            report["arms"].append(row)
            write_json(output / "results.json", report)
            started = time.monotonic()
            verify_inputs(plan, parent, identities)
            arm = output / old["label"]
            arm.mkdir(exist_ok=False)
            (arm / "host_before.txt").write_text(host_activity())
            native = old["native"]
            native_case = case["original_case"]
            if old["action"] == "run_missing":
                native_case = case["fresh_case"]
                native_started = time.monotonic()
                native = run_missing(native_case, parent, env)
                row["fresh_native_wall_s"] = time.monotonic() - native_started
                verify_inputs(plan, parent, identities)
            if old["action"] == "reuse_scored":
                scoring_folder, scoring_id = Path(native_case["output"]), native_case["scoring_id"]
                row.update({k: old["parent_row"][k] for k in ("scores", "metric_files", "conversion", "command")})
            else:
                scoring_folder, scoring_id = Path(case["scoring_output"]), case["scoring_id"]
                scoring_folder.mkdir(exist_ok=False)
                score_started = time.monotonic()
                row.update(score_qfo(repo, Path(native["phylogeny"]), scoring_folder, scoring_id))
                row["scoring_recovery_wall_s"] = time.monotonic() - score_started
            row["readback"] = independent_scores(repo, row, native, scoring_folder, scoring_id)
            identities.extend([*native["checked"], *row["readback"]["checked"]])
            row.update(native=native, native_case=native_case, scoring_output=str(scoring_folder), scoring_id=scoring_id,
                       native_timing_scope="fresh original or missing-arm inference; see native metrics",
                       total_recovery_wall_s=time.monotonic() - started)
            verify_inputs(plan, parent, identities)
            row.update(status="independent_content_readback_verified", inputs_sources_runtime_reverified=True)
            (arm / "host_after.txt").write_text(host_activity())
            write_json(output / "results.json", report)
        report["status"] = "finished_all_arms_independently_verified"
    except Exception as error:
        report.update(status="failed_with_explicit_gaps", error_type=type(error).__name__, error=str(error))
        if report["arms"] and report["arms"][-1]["status"] == "running":
            report["arms"][-1].update(status="failed", error_type=type(error).__name__, error=str(error))
    write_json(output / "results.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    prep = modes.add_parser("prepare")
    for name in ("repo", "source-root", "output", "parent-plan"):
        prep.add_argument("--" + name, type=Path, required=True)
    prep.add_argument("--parent-job", required=True)
    worker = modes.add_parser("run")
    worker.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source_root.resolve(), args.output.resolve(),
                      args.parent_plan.resolve(), args.parent_job) if args.mode == "prepare"
              else run(args.execution_plan.resolve()))
    print(json.dumps({"status": result["status"]}))
    return 0 if result["status"] in {"prepared_waiting_for_terminal_parent", "finished_all_arms_independently_verified"} else 1


if __name__ == "__main__":
    raise SystemExit(main())

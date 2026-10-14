"""Run only the original unrun QfO sensitivity arm after control recovery."""

import argparse
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.report_additional_hillclimb import (
    canonical_recovery_cell, check_plan, shared_native,
)
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores
from benchmark_tools.run_full_search_hillclimb import native_command, runtime, verify_input_inventory
from benchmark_tools.run_membership_hillclimb import DATASETS, argument, score_qfo, write_json
from benchmark_tools.run_qfo_recovered_assessment import require_darwin_work_path
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import sources, timed_stage
from benchmark_tools.validate_profile_runtime import require_profile_runtime

PROTOCOL = "hillclimb_shared_qfo_continuation_20261014.json"
REVISION = "7e7acaed3f3392197c0c498dc60ba63abb726810"


def original_case(plan):
    rows = [case for case in plan["cases"] if case.get("dataset") == "QfO"]
    if (len(rows) != 1 or plan["source_revision"] != REVISION or plan["cpu"] != 16
            or tuple(case["dataset"] for case in plan["cases"]) != DATASETS):
        raise ValueError("Require the original frozen three-dataset sensitivity plan")
    case = rows[0]
    if (case["label"] != "shared_hit1" or case["search_min_total_hits"] != 1
            or case["search_max_candidates_per_query"] != 100 or len(case["inputs"]) != 78):
        raise ValueError("Wrong original QfO sensitivity intervention or inventory")
    command = case["command"]
    expected = native_command(Path(plan["source_root"]), Path(case["fasta"]), Path(case["output"]),
        [argument(command, "--aligner"), argument(command, "--tree-builder")], 100, 16)
    expected.extend(["--search-result-order", "canonical", "--search-min-total-hits", "1"])
    if expected != command:
        raise ValueError("Original native command differs from the prespecified method")
    return case


def relocated_case(plan, output):
    original = original_case(plan)
    case = deepcopy(original)
    folder = output / "native"
    case["output"] = str(folder)
    case["command"][3:5] = [str(folder / "output"), str(folder / "metrics.json")]
    return case


def original_skip(report, plan, identity):
    if (report.get("job_id") != "24242" or report.get("source_revision") != REVISION
            or report.get("execution_plan") != identity
            or report.get("status") != "finished_with_explicit_gaps"
            or tuple(row.get("dataset") for row in report.get("arms", [])) != DATASETS):
        raise ValueError("Require the terminal original sensitivity report with explicit gaps")
    rows = [row for row in report["arms"] if row.get("dataset") == "QfO"]
    if len(rows) != 1:
        raise ValueError("Ambiguous original QfO arm")
    row = rows[0]
    if (row.get("label") != "shared_hit1" or row.get("status") != "not_run_unverified_control"
            or row.get("error_type") != "ValueError"
            or row.get("error") != "Canonical control unavailable or not independently verified"
            or row.get("scores") or any(key in row for key in ("execution", "native", "cpu_backend_probe", "control"))
            or Path(original_case(plan)["output"]).exists()):
        raise ValueError("Original QfO inference was not explicitly unrun; do not duplicate it")
    return row


def terminal_dependencies(accounting):
    rows = [line.split("|") for line in accounting.splitlines() if line]
    if len(rows) != 2 or {row[0] for row in rows} != {"24242", "24246"}:
        raise ValueError("Require unique accounting for both dependency jobs")
    states = {row[0]: row[1:] for row in rows}
    if states["24246"] != ["COMPLETED", "0:0"] or states["24242"] not in (["FAILED", "1:0"], ["COMPLETED", "0:0"]):
        raise ValueError("Require successful recovery and a terminal original sensitivity job")


def short_id(repo, output):
    name = "sh_" + hashlib.sha256(str(output).encode()).hexdigest()[:12]
    work = repo / "qfo_benchmark/w" / ("hc1_" + name)
    require_darwin_work_path(work)
    if work.exists():
        raise FileExistsError("Continuation scoring work exists")
    return name


def prepare(repo, source, output):
    if output.exists() or source == repo:
        raise ValueError("Require fresh continuation output and frozen runner")
    parent_root = repo / "benchmarks/results/shared_hit_pipeline_r1_v1_20261014"
    parent_path = parent_root / "execution_plan.json"
    parent = json.loads(parent_path.read_text())
    case = relocated_case(parent, output)
    verify_input_inventory(case)
    identifier = short_id(repo, output)
    recovery_root = repo / "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014"
    inputs = [record(source / "benchmark_tools" / PROTOCOL), record(parent_path),
              record(recovery_root / "execution_plan.json"),
              *parent["inputs"], *parent["sources"], *parent["runtime"]["files"]]
    source_pins = [*sources(source), record(source / "benchmark_tools/continue_shared_qfo.slurm")]
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    for pin in (*inputs, *source_pins):
        check(pin)
    for pin in source_pins:
        relative = Path(pin["path"]).relative_to(source)
        if relative.suffix != ".so":
            blob = subprocess.check_output(["git", "-C", str(source), "show", revision + ":" + str(relative)])
            if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
                raise ValueError("Continuation runner differs from committed source")
    if runtime() != parent["runtime"]:
        raise ValueError("Original scientific runtime changed")
    require_profile_runtime(Path(parent["source_root"]))
    plan = {"status": "prepared_waiting_for_original_skip_and_recovered_control", "repo": str(repo),
            "source_root": str(source), "source_revision": revision, "output": str(output), "cpu": 16,
            "scientific_source": parent["source_root"], "scientific_revision": REVISION,
            "parent_plan": record(parent_path), "parent_root": str(parent_root),
            "recovery_root": str(recovery_root), "inputs": inputs, "sources": source_pins,
            "runtime": parent["runtime"], "tool_path": parent["tool_path"], "case": case,
            "scoring_id": identifier, "dependency_jobs": ["24242", "24246"],
            "defaults_promoted": False, "goal_completed": False}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", plan)
    return plan


def run(path):
    plan = json.loads(path.read_text())
    repo, source, science, output, parent_root, recovery_root = (Path(plan[k]) for k in (
        "repo", "source_root", "scientific_source", "output", "parent_root", "recovery_root"))
    if Path(__file__).resolve() != source / "benchmark_tools/continue_shared_qfo.py" or (output / "results.json").exists():
        raise ValueError("Require frozen runner and preserve all attempts")
    result = {"status": "validating_original_skip_and_recovered_control", "execution_plan": record(path),
              "job_id": os.environ.get("SLURM_JOB_ID"), "scores": {}, "host_before": host_activity(),
              "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", result)
    try:
        accounting = subprocess.check_output(["sacct", "-j", "24242,24246", "-X", "-n", "-P",
                                             "-o", "JobID,State%30,ExitCode"], text=True)
        terminal_dependencies(accounting)
        check(plan["parent_plan"])
        parent = json.loads(Path(plan["parent_plan"]["path"]).read_text())
        if relocated_case(parent, output) != plan["case"] or short_id(repo, output) != plan["scoring_id"]:
            raise ValueError("Continuation changed original science or scoring namespace")
        parent_path = parent_root / "results.json"
        parent_report = json.loads(parent_path.read_text())
        original_skip(parent_report, parent, plan["parent_plan"])
        recovery_path = recovery_root / "results.json"
        recovery = json.loads(recovery_path.read_text())
        recovery_plan, recovery_pins = check_plan(repo, recovery_root, recovery)
        if record(recovery_root / "execution_plan.json") not in plan["inputs"] or recovery.get("job_id") != "24246":
            raise ValueError("Recovered control differs from the planned dependency")
        control = canonical_recovery_cell(repo, recovery_plan, recovery, recovery_root)
        write_json(output / "original_skip_snapshot.json", parent_report)
        checked = [record(path), record(parent_path), record(recovery_path),
                   record(output / "original_skip_snapshot.json"), *plan["inputs"], *plan["sources"],
                   *plan["runtime"]["files"], *recovery_pins, *control["checked"]]
        for pin in checked:
            check(pin)
        if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
            raise ValueError("Runtime or CUDA visibility changed")
        result.update(control={"report": record(recovery_path), "cell": control},
                      original_skip=record(output / "original_skip_snapshot.json"), dependency_accounting=accounting)
        case = plan["case"]
        verify_input_inventory(case)
        folder = Path(case["output"])
        folder.mkdir(parents=True, exist_ok=False)
        env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(science), "PYTHONHASHSEED": "0",
               "CUDA_VISIBLE_DEVICES": "", "PYTHONFAULTHANDLER": "1", "OMP_NUM_THREADS": "1",
               "NUMBA_NUM_THREADS": "16", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        subprocess.check_output([sys.executable, "-c",
            "from orthohmm.search import engine; from orthohmm.search.prefilter import is_c_prefilter_available; "
            "assert not engine.is_cuda_available() and engine.is_c_available() and is_c_prefilter_available()"],
            cwd=science, env=env, text=True)
        result["status"] = "fresh_original_shared_qfo_running"
        write_json(output / "results.json", result)
        result["execution"] = timed_stage(case["command"], folder, "native", science, env)
        if result["execution"]["returncode"]:
            raise RuntimeError("Original sensitivity inference failed; retain artifacts")
        native = shared_native(case, parent)
        if (native["genes"], native["species"]) != (984137, 78):
            raise ValueError("Wrong corrected QfO gene universe")
        phylo = folder / "output/orthohmm_phylogeny"
        native.update(phylogeny=str(phylo), prediction=str(folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"),
                      checked=[*native["products"].values(), record(folder / "metrics.json"),
                               record(phylo / "provenance_manifest.json"), record(phylo / "reconciliation_summary.json"),
                               *[record(p) for p in sorted((folder / "output/orthohmm_working_res/high_sensitivity_checkpoint").iterdir())]])
        result.update(native=native, status="native_verified_scoring_pending")
        write_json(output / "results.json", result)
        scoring = output / "scoring"
        scoring.mkdir(exist_ok=False)
        result.update(score_qfo(repo, phylo, scoring, plan["scoring_id"]))
        scores, details, score_pins = readback_scores(repo, "QfO", result, native, scoring,
            "ohmm_hillclimb_r1_" + plan["scoring_id"])
        checked.extend([*native["checked"], *score_pins])
        for pin in checked:
            check(pin)
        if runtime() != plan["runtime"]:
            raise ValueError("Runtime changed during continuation")
        result.update(status="original_shared_qfo_native_and_all_six_scores_verified", details=details,
                      metric_vector=endpoint_vector("QfO", scores), checked=checked,
                      deltas={k: v - control["scores"][k] for k, v in endpoint_vector("QfO", scores).items()},
                      scope="Fresh unchanged-science min1/cap100/canonical QfO continuation. Cached recovered min4 control costs are separate; FAS is unseeded.")
    except Exception as error:
        result.update(status="failed_or_unverified_shared_qfo_continuation", error_type=type(error).__name__, error=str(error))
    result["host_after"] = host_activity()
    write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    for name in ("repo", "source", "output"):
        prep.add_argument("--" + name, type=Path, required=True)
    worker = sub.add_parser("run")
    worker.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve()) if args.command == "prepare"
              else run(args.execution_plan.resolve()))
    print(result["status"])
    raise SystemExit(0 if result["status"] in {"prepared_waiting_for_original_skip_and_recovered_control",
        "original_shared_qfo_native_and_all_six_scores_verified"} else 1)

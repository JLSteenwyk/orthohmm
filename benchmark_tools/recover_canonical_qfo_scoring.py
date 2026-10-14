"""Score a completed canonical QfO arm without changing its failed attempt."""

import argparse
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.recover_full_qfo_search import PATH_ERROR
from benchmark_tools.report_additional_hillclimb import check_plan
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
from benchmark_tools.run_canonical_downstream import commands
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import score_qfo, write_json
from benchmark_tools.run_qfo_recovered_assessment import require_darwin_work_path
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import sources

PROTOCOL = "hillclimb_canonical_qfo_scoring_recovery_20261014.json"
REVISION = "402d32323f7166c6135f4c81b849b0abfc97e55c"


def failed_row(report):
    rows = [r for r in report.get("arms", [])
            if (r.get("dataset"), r.get("label")) == ("QfO", "canonical")]
    if (len(rows) != 1 or report.get("source_revision") != REVISION
            or report.get("job_id") != "24240"):
        raise ValueError("Require the unique original canonical QfO arm")
    row = rows[0]
    if (row.get("status") != "failed" or row.get("error_type") != "ValueError"
            or row.get("error") != PATH_ERROR or row.get("scores")
            or row.get("native_execution", {}).get("returncode") != 0
            or row.get("graph_execution", {}).get("returncode") != 0
            or not row.get("outputs") or not row.get("native_summary")):
        raise ValueError("Require successful native inference and only the scoring-path failure")
    return deepcopy(row)


def scoring_id(repo, output):
    name = "cq_" + hashlib.sha256(str(output).encode()).hexdigest()[:12]
    work = repo / "qfo_benchmark/w" / ("hc1_" + name)
    require_darwin_work_path(work)
    if work.exists():
        raise FileExistsError("Recovery scoring work namespace exists")
    return name


def verify_parent(repo, root, report, row):
    plan, pins = check_plan(repo, root, report)
    if failed_row(report) != row or plan["source_revision"] != REVISION or plan["cpu"] != 16:
        raise ValueError("Original failed row or scientific plan changed")
    cases = [case for case in plan["cases"] if case["dataset"] == "QfO"]
    if len(cases) != 1:
        raise ValueError("Require the unique planned QfO case")
    folder = root / "QfO/canonical"
    preflight_path = folder / "preflight.json"
    preflight = json.loads(preflight_path.read_text())
    checkpoint = Path(row["ordering"]["checkpoint"])
    sha = row["ordering"]["manifest"]["sha256"]
    expected = commands(Path(plan["source_root"]), plan, cases[0], folder, checkpoint, sha)
    if (preflight["graph"], preflight["native"]) != expected:
        raise ValueError("Original graph/native command differs from its frozen plan")
    if row["native_execution"]["command"] != expected[1] or row["graph_execution"]["command"] != expected[0]:
        raise ValueError("Executed commands differ from preflight")
    old_id = "canon_" + hashlib.sha256(str(root).encode()).hexdigest()[:16] + "_canonical"
    try:
        require_darwin_work_path(repo / "qfo_benchmark/w" / ("hc1_" + old_id))
    except ValueError as error:
        if str(error) != PATH_ERROR:
            raise
    else:
        raise ValueError("Original scoring-path failure does not reproduce")
    native = replay_native(row, {"command": expected[1]})
    if (native["genes"], native["species"]) != (984137, 78) or native["summary"] != row["native_summary"]:
        raise ValueError("Native universe or summary differs")
    pins.extend([record(preflight_path), row["ordering"]["manifest"],
                 row["candidate"]["membership_constraints"], *native["checked"]])
    for pin in pins:
        check(pin)
    return native, pins


def prepare(repo, source, output, parent_root):
    if output.exists() or source == repo:
        raise ValueError("Require a fresh namespace and frozen recovery source")
    report = json.loads((parent_root / "results.json").read_text())
    row = failed_row(report)
    identifier = scoring_id(repo, output)
    native, pins = verify_parent(repo, parent_root, report, row)
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    source_pins = [*sources(source), record(source / "benchmark_tools/recover_canonical_qfo_scoring.slurm")]
    for pin in source_pins:
        check(pin)
        relative = Path(pin["path"]).relative_to(source)
        if relative.suffix != ".so":
            blob = subprocess.check_output(["git", "-C", str(source), "show", revision + ":" + str(relative)])
            if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
                raise ValueError("Recovery source differs from committed revision")
    output.mkdir(parents=True)
    snapshot = output / "failed_parent_snapshot.json"
    write_json(snapshot, report)
    plan = {"status": "prepared_scoring_only_recovery", "repo": str(repo), "source_root": str(source),
            "output": str(output), "parent_root": str(parent_root), "cpu": 8,
            "source_revision": revision, "scientific_revision": REVISION,
            "parent_snapshot": record(snapshot), "original_failed_row": row,
            "inputs": pins, "sources": source_pins, "runtime": runtime(),
            "protocol": record(source / "benchmark_tools" / PROTOCOL),
            "scoring_id": identifier, "native": native,
            "defaults_promoted": False, "goal_completed": False}
    write_json(output / "execution_plan.json", plan)
    return plan


def run(path):
    plan = json.loads(path.read_text())
    repo, source, output, root = (Path(plan[k]) for k in ("repo", "source_root", "output", "parent_root"))
    if Path(__file__).resolve() != source / "benchmark_tools/recover_canonical_qfo_scoring.py":
        raise ValueError("Execute the frozen recovery runner")
    if (output / "results.json").exists():
        raise FileExistsError("Preserve every previous recovery attempt")
    pins = [record(path), plan["parent_snapshot"], plan["protocol"],
            *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    result = {"status": "scoring_only_recovery_running", "execution_plan": record(path),
              "job_id": os.environ.get("SLURM_JOB_ID"), "scores": {},
              "original_status": "failed", "host_activity_before": host_activity(),
              "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", result)
    started = time.monotonic()
    try:
        for pin in pins:
            check(pin)
        if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
            raise ValueError("Recovery runtime or CUDA visibility changed")
        report = json.loads((root / "results.json").read_text())
        native, original_pins = verify_parent(repo, root, report, plan["original_failed_row"])
        if native != plan["native"] or scoring_id(repo, output) != plan["scoring_id"]:
            raise ValueError("Recovery native inputs or namespace changed")
        result.update(score_qfo(repo, Path(native["phylogeny"]), output, plan["scoring_id"]))
        scores, details, score_pins = readback_scores(repo, "QfO", result, native, output,
            "ohmm_hillclimb_r1_" + plan["scoring_id"])
        checked = [*pins, *original_pins, *score_pins]
        for pin in checked:
            check(pin)
        if runtime() != plan["runtime"] or failed_row(json.loads((root / "results.json").read_text())) != plan["original_failed_row"]:
            raise ValueError("Original completed arm or runtime changed during scoring")
        result.update(status="canonical_qfo_scoring_recovered_independently_verified",
                      details=details, checked=checked, native=native,
                      metric_vector=endpoint_vector("QfO", scores),
                      original_attempt_preserved=True,
                      scope="Scoring-only recovery of unchanged cached downstream inference; not fresh full-pipeline timing. FAS is unseeded.")
    except Exception as error:
        result.update(status="failed_scoring_only_recovery", error_type=type(error).__name__, error=str(error))
    result.update(scoring_recovery_wall_s=time.monotonic() - started, host_activity_after=host_activity())
    write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    for name in ("repo", "source", "output", "parent-root"):
        prep.add_argument("--" + name, type=Path, required=True)
    execution = sub.add_parser("run")
    execution.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve(), args.parent_root.resolve())
              if args.command == "prepare" else run(args.execution_plan.resolve()))
    print(result["status"])
    raise SystemExit(0 if result["status"] in {"prepared_scoring_only_recovery",
                    "canonical_qfo_scoring_recovered_independently_verified"} else 1)

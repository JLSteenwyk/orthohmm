"""Independently read back the candidate score/order intervention."""

import argparse
import json
from pathlib import Path
import shlex
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_installed_ob_clustering import write_json
from benchmark_tools.audit_ob_dependency_replay import candidate_readback
from benchmark_tools.compare_installed_ob_search import partition
from benchmark_tools.audit_installed_orthobench import compare_partitions
from benchmark_tools.probe_ob_candidate_order_scores import LABELS

PLAN_SHA = "cc3bdcef846c68b5d6a2e524a976d68325cfd17e184fdb2934eac618011a9ae8"


def compare_all(groups):
    if tuple(groups) != LABELS:
        raise ValueError("Missing, duplicate or reordered arms")
    return {a+"__vs__"+b: compare_partitions(groups[a], groups[b])
            for i,a in enumerate(LABELS) for b in LABELS[i+1:]}


def validate_summary(report, plan_record, source_record):
    if (report["status"] != "candidate_score_order_factorial_complete"
            or report["plan"] != plan_record or report["source"] != source_record
            or report["accuracy_evaluated"] is not False or report["phylogeny_run"] is not False
            or report["genes"] != 251378 or report["nonself_hits"] != 18235373
            or report["removed_self_hits"] != 251135
            or [r["label"] for r in report["rows"]] != list(LABELS)):
        raise ValueError("Wrong factorial summary scope, inputs or identities")


def audit(repo, directory):
    records = []
    def read(path, sha=None):
        item = record(path)
        if sha is not None and item["sha256"] != sha:
            raise ValueError("Changed pinned evidence")
        records.append(item)
        return json.loads(Path(path).read_text())
    plan_path = directory / "plan.json"
    plan = read(plan_path, PLAN_SHA)
    if directory != Path(plan["output"]) or plan["repo"] != str(repo):
        raise ValueError("Plan directory mismatch")
    if (directory / "failure.json").exists():
        raise ValueError("Failure receipt present")
    submission, started = read(directory / "submission.json"), read(directory / "started.json")
    report = read(directory / "report.json")
    validate_summary(report, record(plan_path), plan["source"])
    scheduler_command = ["sacct", "-j", "22322", "--format=JobID,State,ExitCode,Elapsed", "-n", "-P"]
    scheduler = subprocess.check_output(scheduler_command, text=True, timeout=60)
    root_rows = [line.split("|") for line in scheduler.splitlines() if line.startswith("22322|")]
    if len(root_rows) != 1 or root_rows[0][1:3] != ["COMPLETED", "0:0"]:
        raise ValueError("Scheduler has not confirmed successful completion")
    expected_command = ["env", "-i", "PATH=/usr/bin:/bin", "HOME=/home/bizon", "LANG=C.UTF-8",
                        *[k+"="+v for k,v in plan["env"].items()], plan["runtime"]["executable"], "-I",
                        plan["source"]["path"], "--run", str(plan_path), "--plan-sha256", PLAN_SHA]
    if (submission["job_id"] != "22322" or submission["plan"] != record(plan_path)
            or submission["command"][-2] != "--wrap"
            or shlex.split(submission["command"][-1]) != expected_command
            or started["plan"] != record(plan_path)
            or started["command"] != expected_command[-5:]
            or started["runtime"] != {k:plan["runtime"][k] for k in ("versions", "sources", "native")}):
        raise ValueError("Execution/runtime differs from planned intervention")
    # env -i deliberately drops scheduler variables; submission/accounting bind the job.
    if started["job_id"] is not None:
        raise ValueError("Unexpected inherited scheduler environment")
    records.extend(plan["checked_records"])
    for item in records:
        check(item)
    previous = read(repo / "benchmark_tools/results/ob_dependency_replay_readback_22320.json",
                    "973d1bf10bb0e4e0028e0762cbee36def2c57e4897aa1308c25e231b87f55df4")
    pinned = {r["path"]: r for r in previous["checked_records"]}
    stage_path = repo / "benchmarks/work/ob_dependency_replay_v2_20260926/leiden011/stage_report.json"
    check(pinned[str(stage_path)])
    prior_stage = read(stage_path, pinned[str(stage_path)]["sha256"])
    names_path = Path(plan["checkpoint"]) / "gene_names.txt"
    names = names_path.read_text().splitlines()
    if len(names) != 251378 or len(set(names)) != len(names):
        raise ValueError("Unexpected gene universe")
    universe, groups, rows = set(names), {}, []
    for row in report["rows"]:
        label = row["label"]
        root = directory / label
        if read(root / "result.json") != row:
            raise ValueError("Arm report differs from aggregate")
        work = root / "replay/orthohmm_working_res"
        if (row["output_records"] != [record(p) for p in sorted(work.iterdir())]
                or row["prediction"] != record(work / "phylogeny_candidate_superfamilies.txt")
                or row["candidate_summary"]["parameters"] != prior_stage["candidates"]["parameters"]):
            raise ValueError("Arm output or parameter identity differs")
        evidence = candidate_readback(root, dict(candidates=row["candidate_summary"]), plan["seed"], universe)
        groups[label] = partition(Path(row["prediction"]["path"]), universe)
        records.extend(row["output_records"])
        rows.append(dict(label=label, prediction=row["prediction"], candidate_readback=evidence))
    comparisons = compare_all(groups)
    if comparisons != report["comparisons"]:
        raise ValueError("Reported comparisons do not reproduce")
    historical = next(r for r in previous["checked_records"]
                      if r["sha256"] == "44f201dbc7e4eccf06d9401e5ce1b20fdf6ad523e15ee6007fe1c58f3941842d")
    baselines = {}
    for label, item in (("historical", historical), ("dependency_011_fresh", prior_stage["candidate_partition"])):
        check(item)
        records.append(item)
        baseline = partition(Path(item["path"]), universe)
        baselines[label] = {k:compare_partitions(baseline,v) for k,v in groups.items()}
    records.append(record(directory / "slurm-22322.log"))
    for item in records:
        check(item)
    return dict(status="candidate_score_order_intervention_read_back", source=record(__file__),
                checked_records=records, rows=rows, comparisons=comparisons, baseline_comparisons=baselines,
                scheduler=dict(command=scheduler_command, output=scheduler), accuracy_evaluated=False,
                limitations=["Fixed-seed/fixed-runtime candidate mechanism diagnostic, not final-F1 causality.",
                    "Single execution per arm; no universal determinism or controlled timing claim.",
                    "Merge consistency independently reconstructed; search support not independently recalculated.",
                    "Scheduler linkage uses retained submission/accounting; isolated child has no SLURM_JOB_ID."])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--directory", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    write_json(a.output, audit(a.repo.resolve(), a.directory.resolve()))

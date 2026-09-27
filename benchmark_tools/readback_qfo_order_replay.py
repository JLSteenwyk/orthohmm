"""Independently read back the frozen paired QfO candidate experiment."""

import argparse
from itertools import combinations, zip_longest
import json
from pathlib import Path
import subprocess

from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.run_qfo_order_replay import ARMS, INSTALL_SHA, record, save, validate_plan
from benchmark_tools.verify_ygob_validation import require_completed_job

PLAN_SHA = "bbbb8b04e10d72839dc51de5b1b0484e2bc61a6802a189616cb64f09872c5c6b"
SUBMISSION_SHA = "b8e681bbaf7376709bfc8fcb9754eaf81c02714a07d0a24b2909b87dd2de29de"
JOB = 22328


def partition_difference(left, right):
    left, right = {frozenset(g) for g in left}, {frozenset(g) for g in right}
    changed = left ^ right
    return dict(equal=left == right, left_groups=len(left), right_groups=len(right),
                shared_groups=len(left & right), left_only_groups=len(left - right),
                right_only_groups=len(right - left),
                genes_in_changed_groups=len(set().union(*changed)) if changed else 0)


def semantic_constraints(rows):
    return [(tuple(sorted(r["source_genes"])), tuple(sorted(r["target_genes"]))) for r in rows]


def constraint_difference(left, right):
    left, right = semantic_constraints(left), semantic_constraints(right)
    first = None
    count = 0
    for index, (a, b) in enumerate(zip_longest(left, right)):
        if a != b:
            count += 1
            if first is None:
                first = index
    return dict(equal=left == right, left_constraints=len(left), right_constraints=len(right),
                unequal_positions=count, first_unequal_position=first)


def verify_arm(plan, plan_record, arm, execution, started, complete):
    directory = Path(plan["directory"])
    source = str(Path(plan["repo"]) / "benchmark_tools/run_qfo_order_replay.py")
    expected = ["/usr/bin/time", "-v", "-o", str(directory / (arm + ".time.txt")),
                plan["python"], "-I", source, "worker", "--plan", plan_record["path"],
                "--sha256", PLAN_SHA, "--arm", arm]
    pinned = {r["path"]: r for r in plan["checked_records"]}
    if (execution["arm"] != arm or execution["returncode"] != 0
            or execution["command"] != expected or started["arm"] != arm
            or started["plan"] != plan_record or started["attempts"] != 1
            or started["executable"] != plan["python"] or started["source"] != pinned[source]
            or pinned.get(started["pipeline"]["path"]) != started["pipeline"]
            or complete["arm"] != arm or complete["plan"] != plan_record
            or complete["status"] != "candidate_arm_complete_pending_readback"
            or complete["accuracy_evaluated"] is not False):
        raise ValueError("Arm provenance or scope differs")


def admit(repo, directory):
    accounting = subprocess.check_output(["sacct", "-j", str(JOB), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, JOB)
    if scheduler["AllocCPUS"] != "2":
        raise ValueError("Wrong CPU allocation")
    path = directory / "plan.json"
    plan_record = record(path)
    plan = validate_plan(path, PLAN_SHA)
    if plan["repo"] != str(repo) or plan["directory"] != str(directory):
        raise ValueError("Wrong locations")
    checked = []

    def read(path):
        checked.append(record(path))
        return json.loads(Path(path).read_text())

    frozen = repo / "benchmark_tools/results/qfo_order_replay_submission_22328.json"
    if record(frozen)["sha256"] != SUBMISSION_SHA:
        raise ValueError("Changed frozen submission")
    submission = read(directory / "submission.json")
    if submission != read(frozen) or submission["plan"] != plan_record or submission["job_id"] != str(JOB):
        raise ValueError("Wrong submission")
    if record(submission["protocol"]["path"]) != submission["protocol"]:
        raise ValueError("Changed protocol")
    checked.append(submission["protocol"])
    if record(directory / "install_preparation.json")["sha256"] != INSTALL_SHA:
        raise ValueError("Changed installation preparation")
    if read(directory / "execution_started.json") != dict(plan=plan_record, job_id=str(JOB)):
        raise ValueError("Wrong execution start")
    execution = read(directory / "execution.json")
    if (execution["status"] != "paired_candidates_complete_pending_readback"
            or execution["plan"] != plan_record or execution["job_id"] != str(JOB)
            or execution["accuracy_evaluated"] is not False
            or [r["arm"] for r in execution["arms"]] != list(ARMS)
            or (directory / "execution_failure.json").exists()):
        raise ValueError("Execution incomplete or failed")
    for arm, executed in zip(ARMS, execution["arms"]):
        target = directory / arm
        if (target / "failure.json").exists():
            raise ValueError("Retained arm failure")
        started, complete = read(target / "started.json"), read(target / "complete.json")
        verify_arm(plan, plan_record, arm, executed, started, complete)
        for key, suffix in (("log", ".log"), ("time", ".time.txt")):
            if executed[key] != record(directory / (arm + suffix)):
                raise ValueError("Changed execution log")
            checked.append(executed[key])
        working = target / "orthohmm_working_res"
        outputs = [record(p) for p in sorted(working.iterdir()) if p.is_file()]
        if complete["outputs"] != outputs:
            raise ValueError("Changed arm outputs")
        checked.extend(outputs)
    return dict(scheduler=scheduler, checked_records=checked, plan=plan_record)


def readback(repo, directory, output):
    if output.exists():
        raise FileExistsError(output)
    admitted = admit(repo, directory)
    plan = json.loads((directory / "plan.json").read_text())
    names = (Path(plan["checkpoint"]["path"]).parent / "gene_names.txt").read_text().splitlines()
    universe = set(names)
    if len(names) != len(universe) or len(universe) != plan["expected_genes"]:
        raise ValueError("Invalid universe")
    paths = {"historical": (Path(plan["historical"]["path"]), Path(plan["historical_constraints"]["path"]))}
    for arm in ARMS:
        working = directory / arm / "orthohmm_working_res"
        paths[arm] = (working / "orthohmm_edges_clustered.txt", working / "phylogeny_candidate_merges.json")
    partitions, constraints = {}, {}
    for name, (partition, trace) in paths.items():
        partitions[name] = read_partition(partition, universe)
        constraints[name] = load_membership_constraints(trace, partition)
    contrasts = {}
    for left, right in combinations(paths, 2):
        contrasts[left + "_vs_" + right] = dict(
            partition=partition_difference(partitions[left], partitions[right]),
            constraints=constraint_difference(constraints[left], constraints[right]))
    if admit(repo, directory) != admitted:
        raise ValueError("Admission changed during readback")
    result = dict(status="paired_candidate_readback_complete", job_id=JOB,
        admission=admitted, genes=len(universe), contrasts=contrasts,
        sources=[record(Path(__file__)), *[record(repo / "benchmark_tools" / name) for name in
            ("audit_historical_profile_ablation.py", "replay_phylogeny.py", "verify_ygob_validation.py")]],
        accuracy_evaluated=False, historical_scores_replaced=False, publication_ready=False,
        limitations=["Cached candidate stage only; not end-to-end equivalence or accuracy evidence",
                     "Shared-host runtime is not a controlled efficiency comparison"])
    save(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(readback(args.repo.resolve(), args.directory.resolve(), args.output.resolve()), indent=2))

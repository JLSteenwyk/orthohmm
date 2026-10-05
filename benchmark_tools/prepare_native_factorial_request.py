"""Prepare one held request for the next frozen identity; never submit or release."""

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.probe_host_counters import snapshot
from benchmark_tools.run_native_factorial_cost import (
    MEMORY, SCRIPT, available_memory, read, reviewed_history, validate_plan, validate_request,
)

PREPARATION = ROOT / "benchmark_tools/results/native_factorial_receipt_amendment_20261004/preparation.json"
PREPARATION_SHA256 = "30bfac83243891ce494f893e206779c0e4315ef32b9f18ea63d4b7f05dacb314"


def pinned(path, sha256):
    ref = record(path)
    if ref["sha256"] != sha256:
        raise ValueError("Explicit evidence checksum differs: " + str(path))
    return ref


def held_job(raw, job):
    pairs = re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", raw)
    fields = dict(pairs)
    expected = dict(JobId=str(job), JobName="orthohmm_factorial_cost", JobState="PENDING",
        Reason="JobHeldUser", Partition="gpu", ReqNodeList="bizon", NumCPUs="64", NumTasks="1",
        MinMemoryNode="128G", Requeue="0", Restarts="0", TimeLimit="1-02:00:00",
        Command=str(SCRIPT), WorkDir=str(ROOT))
    expected["CPUs/Task"] = "64"
    if len(fields) != len(pairs) or any(fields.get(k) != v for k, v in expected.items()):
        raise ValueError("Held job identity or resource envelope differs")
    return fields


def prior_score(ref, prior_ref, prior):
    score = read(ref)
    if (score.get("schema") != "native_factorial_orthobench_score_v1"
            or score.get("status") != "terminal_native_orthobench_scored"
            or score.get("terminal_review") != prior_ref
            or any(score.get(k) != prior.get(k) for k in ("index", "cell", "job_id", "plan"))
            or score.get("native_outputs_validated") is not True
            or score.get("accuracy_evaluated") is not True):
        raise ValueError("Previous successful OrthoBench output lacks its bound separate score")


def prepare(job, history, score_ref, output):
    output = Path(output)
    if (not output.is_absolute() or output.resolve() != output or output.is_symlink()
            or not output.parent.is_dir() or output.exists()):
        raise ValueError("Require a fresh direct request path in an existing directory")
    preparation_ref = pinned(PREPARATION, PREPARATION_SHA256)
    preparation = read(preparation_ref)
    plan_ref, policy_ref = preparation["plan"], preparation["policy"]
    plan, policy = read(plan_ref), read(policy_ref)
    runs = validate_plan(plan)
    index = len(history)
    if not 3 <= index < len(runs) or policy["plan_sha256"] != plan_ref["sha256"]:
        raise ValueError("Require the remaining reviewed-prefix identity and unchanged policy")
    checked_history = reviewed_history(dict(index=index, history=history), plan_ref, plan)
    priors = [read(ref) for ref in history]
    if type(job) is not int or job <= max(prior["job_id"] for prior in priors):
        raise ValueError("Require a new held job, never a previous attempt")
    prior = priors[-1]
    evidence = [preparation_ref, record(__file__)]
    if prior["dataset"] == "orthobench" and prior["status"] == "native_success":
        if score_ref is None:
            raise ValueError("Require the preceding OrthoBench score before preparing its successor")
        prior_score(score_ref, history[-1], prior)
        evidence.append(score_ref)
    elif score_ref is not None:
        raise ValueError("Do not substitute an OrthoBench score for a failed attempt or QfO endpoint")
    for ref in [*plan["helper_sources"], *plan["evidence"], *policy["evidence"], *evidence]:
        check(ref)
    for path in (Path(runs[index]["output_root"]),
                 Path(plan["panel_root"]) / "sessions" / f"run_{index:02d}"):
        if path.exists() or path.is_symlink():
            raise ValueError("Next identity already has output; no implicit retry or resume")
    raw_memory = Path("/proc/meminfo").read_text()
    capacity = available_memory(raw_memory)
    if capacity < MEMORY:
        raise ValueError("Unsafe current available RAM")
    command = ["scontrol", "show", "job", str(job), "--oneliner"]
    held = subprocess.run(command, capture_output=True, text=True, check=True, timeout=5)
    held_job(held.stdout, job)
    request = dict(schema="native_factorial_cost_request_v1", execution_authorized=True,
        job_id=job, index=index, plan=plan_ref, policy=policy_ref, history=history,
        scheduler_command=str(SCRIPT), allocation_cwd=str(ROOT),
        source_commit=subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        prepared_unix_ns=time.time_ns(), evidence=evidence,
        authorization_basis="Active publication goal; next different frozen identity after explicit reviewed prefix; shared-host authorization",
        prior_history_scheduler_checks=checked_history,
        held_scheduler=dict(command=command, stdout=held.stdout, stderr=held.stderr),
        capacity_precheck=dict(raw_meminfo=raw_memory, available_memory_bytes=capacity,
            host_counters=snapshot(), capacity_guaranteed_through_run=False,
            background_cpu_used_for_eligibility=False),
        automatic_retry=False, timing_success_established=False)
    validate_request(request, plan_ref, job)
    for ref in [preparation_ref, plan_ref, policy_ref, *history, *evidence]:
        check(ref)
    save(output, request)
    return record(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--history", nargs=2, action="append", required=True,
                        metavar=("REVIEW", "SHA256"))
    parser.add_argument("--prior-score", nargs=2, metavar=("SCORE", "SHA256"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    history = [pinned(Path(path), sha256) for path, sha256 in args.history]
    score_ref = pinned(Path(args.prior_score[0]), args.prior_score[1]) if args.prior_score else None
    print(json.dumps(prepare(args.job, history, score_ref, args.output), sort_keys=True))


if __name__ == "__main__":
    main()

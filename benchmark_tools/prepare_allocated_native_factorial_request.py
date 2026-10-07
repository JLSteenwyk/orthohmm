"""Bind one held next-identity request; no submission, release or retry."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import time

from benchmark_tools.native_factorial_allocated_execution import (
    ROOT, SCRIPT, MEMORY, amendment, validate_request, reviewed_history)
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_native_factorial_request import PREPARATION, PREPARATION_SHA256, pinned
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.probe_host_counters import snapshot
from benchmark_tools.run_native_factorial_cost import read, available_memory


def held_job(raw, job):
    pairs = re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", raw)
    fields = dict(pairs)
    expected = dict(JobId=str(job), JobName="orthohmm_allocated_factorial", JobState="PENDING",
        Reason="JobHeldUser", Partition="gpu", ReqNodeList="bizon", NumCPUs="64", NumTasks="1",
        MinMemoryNode="128G", Requeue="0", Restarts="0", TimeLimit="1-02:00:00",
        Command=str(SCRIPT), WorkDir=str(ROOT))
    expected["CPUs/Task"] = "64"
    if (type(job) is not int or job <= 0 or len(fields) != len(pairs)
            or any(fields.get(k) != v for k, v in expected.items())
            or fields.get("NumNodes") not in {"1", "1-1"}
            or not fields.get("UserId", "").endswith("(" + str(os.getuid()) + ")")
            or any(k in fields for k in ("ArrayJobId", "ArrayTaskId", "HetJobId", "HetJobOffset"))):
        raise ValueError("Held allocated-native ownership, identity or envelope differs")
    return fields


def prepare(job, amendment_ref, history, output):
    output = Path(output)
    if (not output.is_absolute() or output.resolve() != output or output.is_symlink()
            or not output.parent.is_dir() or output.exists()):
        raise ValueError("Require a fresh direct request path")
    execution, plan = amendment(amendment_ref)
    preparation_ref = pinned(PREPARATION, PREPARATION_SHA256)
    preparation = read(preparation_ref)
    plan_ref, policy_ref = execution["historical_plan"], preparation["policy"]
    if preparation["plan"] != plan_ref or read(policy_ref)["plan_sha256"] != plan_ref["sha256"]:
        raise ValueError("Shared-host preparation/policy differs")
    index = len(history)
    request = dict(schema="allocated_native_factorial_request_v1", execution_authorized=True,
        job_id=job, index=index, plan=plan_ref, amendment=amendment_ref, policy=policy_ref,
        history=history, scheduler_command=str(SCRIPT), allocation_cwd=str(ROOT), automatic_retry=False)
    validate_request(request, amendment_ref, execution, job)
    checked_history = reviewed_history(request, amendment_ref, execution, plan)
    if job <= max(read(ref)["job_id"] for ref in history):
        raise ValueError("Require a new held job, not a previous attempt")
    run = plan["runs"][index]
    for path in (Path(run["output_root"]), Path(plan["panel_root"]) / "sessions" / f"run_{index:02d}"):
        if path.exists() or path.is_symlink():
            raise ValueError("Next identity already has output; no implicit retry/resume")
    meminfo = Path("/proc/meminfo").read_text()
    capacity = available_memory(meminfo)
    if capacity < MEMORY:
        raise ValueError("Unsafe available RAM")
    query = ["scontrol", "show", "job", str(job), "--oneliner"]
    result = subprocess.run(query, check=True, capture_output=True, text=True, timeout=5)
    held_job(result.stdout, job)
    evidence = [preparation_ref, amendment_ref, record(__file__)]
    request.update(source_commit=subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        prepared_unix_ns=time.time_ns(), evidence=evidence, prior_history_scheduler_checks=checked_history,
        authorization_basis="Active goal; new unrun frozen identity after reviewed prefix and complete placement amendment",
        held_scheduler=dict(command=query, stdout=result.stdout, stderr=result.stderr),
        capacity_precheck=dict(raw_meminfo=meminfo, available_memory_bytes=capacity, host_counters=snapshot(),
            capacity_guaranteed_through_run=False, background_cpu_used_for_eligibility=False),
        timing_success_established=False)
    for ref in [plan_ref, policy_ref, *history, *evidence]:
        check(ref)
    save(output, request)
    return record(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--amendment", type=Path, required=True)
    parser.add_argument("--amendment-sha256", required=True)
    parser.add_argument("--history", nargs=2, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ref = pinned(args.amendment, args.amendment_sha256)
    history = [pinned(Path(path), digest) for path, digest in args.history]
    print(json.dumps(prepare(args.job, ref, history, args.output.absolute()), sort_keys=True))

"""Bind the already held 22428; never submit/release or retry inference."""

import json
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.probe_host_counters import snapshot
from benchmark_tools.run_native_factorial_cost import MEMORY, SCRIPT, available_memory, read, reviewed_history, validate_request


def main():
    folder = ROOT / "benchmark_tools/results/native_factorial_receipt_amendment_20261004"
    preparation_ref = record(folder / "preparation.json")
    if preparation_ref["sha256"] != "30bfac83243891ce494f893e206779c0e4315ef32b9f18ea63d4b7f05dacb314":
        raise ValueError("Prepared amendment identity changed")
    preparation = read(preparation_ref)
    plan_ref, policy_ref = preparation["plan"], preparation["policy"]
    plan, policy = read(plan_ref), read(policy_ref)
    if policy["plan_sha256"] != plan_ref["sha256"] or plan["runs"][1]["cell"] != "p0_c0_r1":
        raise ValueError("Wrong plan/policy/next identity")
    for ref in [*plan["helper_sources"], *plan["evidence"], *policy["evidence"]]:
        check(ref)
    history = [plan["history_adoption"]["prior_review"]]
    checked = reviewed_history(dict(index=1, history=history), plan_ref, plan)
    raw_memory = Path("/proc/meminfo").read_text()
    capacity = available_memory(raw_memory)
    if capacity < MEMORY:
        raise ValueError("Unsafe current available RAM")
    command = ["scontrol", "show", "job", "22428", "--oneliner"]
    held = subprocess.run(command, capture_output=True, text=True, check=True, timeout=5)
    pairs = re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", held.stdout)
    fields = dict(pairs)
    expected = dict(JobId="22428", JobName="orthohmm_factorial_cost", JobState="PENDING", Reason="JobHeldUser",
        Partition="gpu", ReqNodeList="bizon", NumCPUs="64", NumTasks="1", MinMemoryNode="128G",
        Requeue="0", Restarts="0", TimeLimit="1-02:00:00", Command=str(SCRIPT), WorkDir=str(ROOT))
    expected["CPUs/Task"] = "64"
    if len(fields) != len(pairs) or any(fields.get(k) != v for k, v in expected.items()):
        raise ValueError("Held scheduler identity or resource envelope differs")
    source = record(__file__)
    request = dict(schema="native_factorial_cost_request_v1", execution_authorized=True, job_id=22428,
        index=1, plan=plan_ref, policy=policy_ref, history=history, scheduler_command=str(SCRIPT),
        allocation_cwd=str(ROOT), source_commit=subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        prepared_unix_ns=time.time_ns(), evidence=[preparation_ref, source,
            record(ROOT / "benchmark_tools/results/native_factorial_receipt_amendment_tests_20261004.xml")],
        authorization_basis="Active publication goal; next different frozen identity after explicit reviewed history adoption; shared-host authorization",
        prior_history_scheduler_checks=checked, held_scheduler=dict(command=command, stdout=held.stdout, stderr=held.stderr),
        capacity_precheck=dict(raw_meminfo=raw_memory, available_memory_bytes=capacity, host_counters=snapshot(),
            capacity_guaranteed_through_run=False, background_cpu_used_for_eligibility=False),
        automatic_retry=False, timing_success_established=False)
    validate_request(request, plan_ref, 22428)
    output = Path(__file__).with_name("request_01_receipt_amended.json")
    save(output, request)
    print(json.dumps(record(output), sort_keys=True))


if __name__ == "__main__":
    main()

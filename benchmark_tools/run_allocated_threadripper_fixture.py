"""One prepared scheduler-placement fixture; never native inference or retry."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

from benchmark_tools.measure_allocated_threadripper_scaling import measure
from benchmark_tools.replay_allocated_threadripper_scaling import replay
from benchmark_tools.derive_threadripper_resources import endpoints
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_native_factorial_cost import validate_plan, read, available_memory
from benchmark_tools.verify_threadripper_controller import validate as validate_controller
from benchmark_tools.verify_lineage_native_provenance import same
from benchmark_tools.measure_allocated_threadripper_scaling import sources as collector_sources


ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "benchmark_tools/results/allocated_threadripper_fixture_20261006.sh"
PYTHON = ROOT / "benchmarks/work/native_factorial_review_py310_20261004/bin/python"
PLAN = ROOT / "benchmark_tools/results/native_factorial_receipt_amendment_20261004/plan.json"
PLAN_SHA = "6c87babcbb5581830e0b9e7b9bf9aaba30a85bde1c4ab465e561017e67e9c89b"
PROGRAM = """import concurrent.futures, json, os, time
def work(index):
    affinity = sorted(os.sched_getaffinity(0))
    data = bytearray(4*1024**2)
    start = time.monotonic()
    value = index + 1
    while time.monotonic() - start < 6:
        for _ in range(1000):
            value = (value * 1664525 + 1013904223) % 4294967296
    return dict(index=index, pid=os.getpid(), affinity=affinity, checksum=value,
                allocated_bytes=len(data), cgroup=open('/proc/self/cgroup').read())
affinity = sorted(os.sched_getaffinity(0))
with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
    workers = list(pool.map(work, range(4)))
assert len(affinity) == 32 and all(row['affinity'] == affinity for row in workers)
print(json.dumps(dict(status='engineering_work_completed', affinity=affinity, workers=workers)))
"""


def command():
    return [str(PYTHON), "-I", "-S", "-B", "-c", PROGRAM]


def sources():
    paths = [ROOT / "benchmark_tools" / name for name in (
        "run_allocated_threadripper_fixture.py", "replay_allocated_threadripper_scaling.py",
        "replay_threadripper_scaling.py", "derive_threadripper_resources.py",
        "summarize_fixture_memory_scopes.py")]
    return [*collector_sources(), *(record(path) for path in paths), record(SCRIPT), record(PYTHON)]


def budget(raw, job):
    result = validate_controller(raw, job, "running", command=str(SCRIPT), cwd=str(ROOT),
                                 time_limit="00:20:00", allocation_mode="shared")
    elapsed = result["fields"].get("RunTime", "")
    match = re.fullmatch(r"00:([0-5]\d):([0-5]\d)", elapsed)
    if match is None or 1200 - (int(match[1])*60 + int(match[2])) < 300:
        raise ValueError("Insufficient fixture release/reporting allowance")
    return result


def execute(prepared_path, prepared_sha):
    prepared_ref = record(prepared_path)
    if prepared_ref["sha256"] != prepared_sha:
        raise ValueError("Prepared fixture digest differs")
    prepared = read(prepared_ref)
    if (prepared.get("schema") != "allocated_threadripper_fixture_prepared_v1"
            or prepared.get("root") != str(ROOT) or prepared.get("command") != command()
            or prepared.get("native_inference_authorized") is not False
            or prepared.get("automatic_retry") is not False
            or not same(prepared.get("new_sources"), sources())
            or not same(prepared.get("historical_plan"), record(PLAN))
            or prepared["historical_plan"]["sha256"] != PLAN_SHA):
        raise ValueError("Require prepared engineering-only fixture")
    for item in prepared["new_sources"]:
        check(item)
    plan = read(prepared["historical_plan"])
    validate_plan(plan)
    for item in plan["helper_sources"]:
        check(item)
    job = int(os.environ["SLURM_JOB_ID"])
    output = Path(prepared["output_directory"])
    if (not output.is_absolute() or output.resolve() != output
            or not output.is_relative_to(ROOT / "benchmarks/work")):
        raise ValueError("Require direct fresh engineering output under work root")
    output.mkdir(exist_ok=False)
    def guard(directory):
        started = time.monotonic()
        result = subprocess.run(["scontrol", "show", "job", str(job), "--oneliner"],
            check=True, capture_output=True, text=True, timeout=5)
        verified = budget(result.stdout, job)
        if (verified["fields"].get("Comment") != "ohmm-allocated-fixture:" + prepared_sha
                or verified["fields"].get("JobName") != "ohmm_allocated_core_fixture"):
            raise ValueError("Fixture controller differs from held prepared submission")
        meminfo = Path("/proc/meminfo").read_text()
        memory = available_memory(meminfo)
        if memory < 128*1024**3 or time.monotonic() - started > 5:
            raise ValueError("Unsafe memory capacity or stale fixture release observation")
        save(directory / "fixture_release_check.json", dict(controller=verified,
            raw=result.stdout, meminfo=meminfo, memory_available_bytes=memory,
            job_id=job, scientific_execution_authorized=False))
    measured = measure(command(), output / "measurement", job, 32, 128*1024**3,
        85800, 1., release_guard=guard)
    replayed = replay(output / "measurement", job, command())
    save(output / "replay.json", replayed)
    derived = endpoints(replayed, measured["native"], job)
    save(output / "resources.json", derived)
    native = json.loads((output / "measurement/native.log").read_text())
    allowed = replayed["native_cpu_ids"]
    if (native.get("status") != "engineering_work_completed" or native.get("affinity") != allowed
            or len(native.get("workers", [])) != 4
            or any(row["affinity"] != allowed or row["cgroup"] != measured["placement"]["cgroup"]
                   for row in native["workers"])
            or replayed["native_outcome"] != "exited_zero"
            or any(status != "observed_within_affinity" for status in replayed["affinity_observation_statuses"])):
        raise ValueError("Fixture workload/affinity does not reproduce")
    for item in [prepared_ref, *prepared["new_sources"], *plan["helper_sources"], *replayed["evidence"]]:
        check(item)
    result = dict(schema="allocated_threadripper_fixture_result_v1", job_id=job,
        status="engineering_fixture_replayed_pending_terminal_check", prepared=prepared_ref,
        source=record(__file__), replay=record(output / "replay.json"), resources=record(output / "resources.json"),
        native_cpu_ids=allowed, native_log=record(output / "measurement/native.log"),
        native_inference_performed=False, accuracy_admitted=False, production_execution_authorized=False,
        scientific_timings_admitted=False, publication_ready=False,
        limitations=["Short allocation-aware engineering fixture, not production inference or scaling.",
            "Scheduler terminal outcome and external submission/release binding still require independent checks.",
            "Native controller, reviewer and conversion/admission integration remain mandatory before production."])
    save(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--prepared-sha256", required=True)
    args = parser.parse_args()
    print(execute(args.prepared.absolute(), args.prepared_sha256)["status"])

"""Record fresh22429 native startup only; do not admit scores or final costs."""

import json
from pathlib import Path
import subprocess
import sys
import time

import psutil

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.review_native_factorial_attempt import native_command
from benchmark_tools.run_native_factorial_cost import SCRIPT, read, validate_plan, validate_request, MEMORY
from benchmark_tools.verify_threadripper_controller import validate


def main():
    source = record(__file__)
    request_ref = record(Path(__file__).with_name("request_02_receipt_amended.json"))
    if request_ref["sha256"] != "ed7f7850df24c6f90a18c5bda4fbcb4c1f5ad34d5aa2538f879c6decde6312cc":
        raise ValueError("Request checksum differs")
    request = read(request_ref)
    validate_request(request, request["plan"], 22429)
    plan = read(request["plan"])
    run = validate_plan(plan)[2]
    root = Path(run["output_root"])
    receipt = json.loads((root / "native_execution.json").read_text())
    if receipt.get("status") != "native_factorial_running" or receipt.get("index") != 2 or receipt.get("plan") != request["plan"]:
        raise ValueError("Require this native running receipt")
    proc = psutil.Process(receipt["placement"]["pid"])
    command = native_command(request["plan"], run, read(plan["baseline"]))
    affinity = proc.cpu_affinity()
    if proc.cmdline() != command or affinity != list(range(32)) or receipt["placement"]["affinity"] != affinity:
        raise ValueError("Actual native process/affinity differs")
    if len({(r["package"], r["core"]) for r in receipt["placement"]["topology"]}) != 32:
        raise ValueError("Native CPUs are not distinct physical cores")
    if Path(f"/proc/{proc.pid}/cgroup").read_text() != receipt["placement"]["cgroup"]:
        raise ValueError("Actual native cgroup differs")
    limits = [Path(r["path"]) / "memory.max" for r in receipt["placement"]["ancestors"]]
    if str(MEMORY) not in {p.read_text().strip() for p in limits}:
        raise ValueError("Native ancestor has no128GiB cap")
    scheduler_command = ["scontrol", "show", "job", "22429", "--oneliner"]
    done = subprocess.run(scheduler_command, capture_output=True, text=True, check=True, timeout=5)
    validated = validate(done.stdout, 22429, "running", command=str(SCRIPT), cwd=str(ROOT),
        time_limit="1-02:00:00", allocation_mode="shared")
    if validated["fields"]["Comment"] != request_ref["sha256"]:
        raise ValueError("Current scheduler request comment differs")
    paths = [Path(plan["panel_root"]) / "sessions/run_02/started.json",
        Path(plan["panel_root"]) / "sessions/run_02/lookup_checks/checked_01.json",
        root / "preparation.json", *[root / "measurement" / name for name in (
            "environment_preflight.json", "release_budget.json", "launch_environment_observation.json", "command.json")]]
    refs = [record(p) for p in paths]
    preflight = read(refs[3])
    if preflight["decision"] != "passed" or preflight["job_id"] != 22429 or preflight["index"] != 2:
        raise ValueError("Native environmental preflight differs")
    launch = read(refs[5])
    if launch["available_memory_bytes"] < MEMORY or not proc.is_running() or proc.cmdline() != command:
        raise ValueError("Unsafe recorded launch capacity or native process no longer live")
    for ref in [source, request_ref, request["plan"], *refs]:
        check(ref)
    result = dict(schema="native_factorial_live_start_observation_v1", source=source,
        observed_unix_ns=time.time_ns(), job_id=22429, index=2, cell=run["cell"], request=request_ref, plan=request["plan"],
        scheduler=dict(command=scheduler_command, stdout=done.stdout, stderr=done.stderr, validated=validated),
        native_process=dict(pid=proc.pid, create_time=proc.create_time(), command=command, affinity=affinity,
            cpu_times_at_observation=list(proc.cpu_times())), running_receipt_snapshot=receipt, evidence=refs,
        launch_available_memory_bytes=launch["available_memory_bytes"], terminal_reviewed=False,
        accuracy_evaluated=False, scientific_timings_admitted=False, next_identity_authorized=False,
        automatic_retry=False, uncontended_timing=False,
        limitations=["Live startup observation only, not terminal accuracy/resources or continuous integrity.",
            "Shared-host contention has unknown potentially tool-dependent timing effects."])
    output = ROOT / "benchmark_tools/results/native_factorial_native_start_22429.json"
    save(output, result)
    print(json.dumps(record(output), sort_keys=True))


if __name__ == "__main__":
    main()

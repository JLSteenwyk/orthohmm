"""Long-run root-context collection; not an execution authorization or admission.

Keep historical 900-second collectors unchanged. A separately frozen scaling
recipe must select this collector and its matching long-run replay explicitly.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.measure_native_interval_step import run_command, step_memory
from benchmark_tools.measure_native_hierarchy_step import interval_point
from benchmark_tools.measure_native_lineage_step import evaluate as evaluate_lineage
from benchmark_tools.measure_native_root_context import read_point, evaluate, lineage_identity
from benchmark_tools.probe_host_counters import snapshot
from benchmark_tools.probe_dgx_step_separation import save, wait_file

TIMEOUT = 85800


def validate(command, cpus, timeout_s, interval_s):
    if (not isinstance(command, list) or not command
            or not all(isinstance(v, str) and v for v in command)
            or not Path(command[0]).is_absolute()):
        raise ValueError("Require nonempty command with absolute executable")
    if (type(cpus) is not int or cpus != 20 or type(timeout_s) is not int or timeout_s != TIMEOUT
            or type(interval_s) not in (int, float) or interval_s != 1.):
        raise ValueError("Require frozen 20CPU/85800s/1s scaling settings")


def worker(directory):
    plan = json.loads((directory / "command.json").read_text())
    validate(plan["command"], plan["cpus"], plan["timeout_s"], plan["interval_s"])
    save(directory / "ready.json", dict(pid=os.getpid(), cgroup=Path("/proc/self/cgroup").read_text()))
    gate = wait_file(directory / "go.json")
    if not isinstance(gate, dict) or set(gate) != {"go"} or gate["go"] is not True:
        save(directory / "aborted_before_native.json", dict(status="observer_did_not_release_native"))
        return
    before = snapshot()
    start = time.monotonic_ns()
    with (directory / "native.log").open("x") as log:
        code, timed_out = run_command(plan["command"], log, TIMEOUT)
    finish = time.monotonic_ns()
    after = snapshot()
    save(directory / "done.json", dict(exit_code=code, timed_out=timed_out,
        started_ns=start, finished_ns=finish, snapshots=[before, after]))
    wait_file(directory / "release.json")


def measure(command, directory, job_id, cpus, memory_bytes, timeout_s, interval_s,
            monitor_host=True, host_interval_s=30.):
    validate(command, cpus, timeout_s, interval_s)
    if (type(memory_bytes) is not int or memory_bytes != 96*1024**3
            or monitor_host is not True or type(host_interval_s) not in (int, float) or host_interval_s != 30.):
        raise ValueError("Require fixed memory and observer settings")
    if (type(job_id) is not int or job_id <= 0 or os.environ.get("SLURM_JOB_ID") != str(job_id)
            or os.environ.get("SLURM_CPUS_PER_TASK") != "20"
            or os.environ.get("SLURM_MEM_PER_NODE") != "98304" or os.uname().nodename != "spark-7ff0"):
        raise ValueError("Require matching DGX20CPU/96GiB allocation")
    directory = Path(directory).absolute()
    directory.mkdir(exist_ok=False)
    save(directory / "command.json", dict(command=command, cpus=cpus, timeout_s=timeout_s, interval_s=interval_s))
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
                sys.executable, "-B", str(Path(__file__).resolve()), "--worker", str(directory)]
    with (directory / "step.log").open("x") as log:
        process = subprocess.Popen(launched, stdout=log, stderr=subprocess.STDOUT)
        try:
            ready = wait_file(directory / "ready.json")
            points = [read_point(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json")]
            save(directory / "point_000000.json", points[0])
            save(directory / "go.json", {"go": True})
            start = time.monotonic()
            index = 0
            while True:
                index += 1
                time.sleep(max(0, start + index * interval_s - time.monotonic()))
                completed = (directory / "done.json").exists()
                if process.poll() is not None:
                    raise RuntimeError("Worker exited before final observation")
                points.append(read_point(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json"))
                save(directory / f"point_{index:06d}.json", points[-1])
                if completed:
                    break
                if time.monotonic() - start > TIMEOUT + 30:
                    raise TimeoutError("Native command exceeded timeout and cleanup allowance")
            done = json.loads((directory / "done.json").read_text())
            memory = step_memory(interval_point(points[-1], job_id))
            save(directory / "step_memory.json", memory)
            save(directory / "release.json", {"release": True})
            if process.wait(timeout=45) != 0:
                raise RuntimeError("Native step wrapper failed")
            measured = dict(status="command_exited_zero" if done["exit_code"] == 0 else "command_failed",
                native=done, native_wall_s=(done["finished_ns"]-done["started_ns"])/1e9,
                job_id=job_id, launched=launched, points=points, step_memory=memory,
                screening=evaluate_lineage(points, done, job_id), scientific_timings_admitted=False,
                controlled_workload_verified=False, publication_ready=False,
                limitations=["Long-run collector, not controlled comparative timing admission.",
                    "All CPU flags and non-atomic read windows retained; no overhead subtraction.",
                    "Host reads are counter snapshots, not a full process/GPU/device-I/O inventory.",
                    "Point retention and final evaluation consume observer resources; long-run overhead remains unvalidated."])
            save(directory / "lineage_report.json", measured)
            context = dict(status="native_root_context_measured", job_id=job_id,
                native_wall_s=measured["native_wall_s"], context=evaluate(points, job_id),
                lineage_report=lineage_identity(directory), scientific_timings_admitted=False,
                environmental_validity_established=False)
            save(directory / "root_context_report.json", context)
            return measured
        finally:
            # A failed initial observation must not release an unobserved native command.
            if not (directory / "go.json").exists():
                save(directory / "go.json", {"abort": True})
            if not (directory / "release.json").exists():
                save(directory / "release.json", {"release": True})
            if process.poll() is None:
                process.wait(timeout=TIMEOUT + 90)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, required=True)
    worker(parser.parse_args().worker.resolve())

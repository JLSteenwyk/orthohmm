"""Boundary-only native-point control; no launch CLI or timing admission."""

import json
import os
from pathlib import Path
import subprocess
import sys
import time

from benchmark_tools import measure_threadripper_scaling as periodic
from benchmark_tools.command_host_monitor import HostMonitor
from benchmark_tools.disk_observation_sequence import DiskObservations
from benchmark_tools.native_completion import completion_evidence
from benchmark_tools.observe_threadripper_process_identity import enriched_snapshot
from benchmark_tools.periodic_host_observer import PeriodicHostObserver
from benchmark_tools.probe_dgx_step_separation import save, wait_file
from benchmark_tools.report_finalization import observe as observe_reporting
from benchmark_tools.replay_threadripper_scaling import native_outcome
from benchmark_tools.slurm_resource_snapshot import scoped_path

SCHEMA = "threadripper_boundary_control_v1"
POLICY = dict(native_points=2, periodic_native_sampling=False,
              completion_poll_interval_s=1., common_host_interval_s=30.)


def measure(command, directory, job_id, cpus, memory_bytes, timeout_s, interval_s,
            monitor_host=True, host_interval_s=30., *, release_guard=None):
    periodic.validate(command, cpus, timeout_s, interval_s)
    if release_guard is not None and not callable(release_guard):
        raise ValueError("Release guard must be callable")
    if (type(memory_bytes) is not int or memory_bytes != 128 * 1024**3
            or monitor_host is not True or type(host_interval_s) not in (int, float)
            or host_interval_s != 30.):
        raise ValueError("Require fixed memory and common host observer settings")
    if (type(job_id) is not int or job_id <= 0 or os.environ.get("SLURM_JOB_ID") != str(job_id)
            or os.environ.get("SLURM_CPUS_PER_TASK") != "64"
            or os.environ.get("SLURM_MEM_PER_NODE") != "131072" or os.uname().nodename != "bizon"):
        raise ValueError("Require matching Threadripper64-slot/128GiB allocation")
    directory = Path(directory).absolute()
    directory.mkdir(exist_ok=False)
    save(directory / "command.json", dict(command=command, cpus=cpus,
        timeout_s=timeout_s, interval_s=interval_s))
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1",
        "--cpus-per-task=64", "--cpu-bind=mask_cpu:0xffffffff", sys.executable,
        "-B", str(Path(periodic.__file__).resolve()), "--worker", str(directory)]
    with (directory / "step.log").open("x") as log, (directory / "host_processes.jsonl").open("x") as host_log:
        process = subprocess.Popen(launched, stdout=log, stderr=subprocess.STDOUT)
        host_observer = None
        try:
            ready = wait_file(directory / "ready.json")
            scope = scoped_path(ready["cgroup"], job_id)
            job_scope = next(p for p in scope.parents if p.name == f"job_{job_id}")
            host = HostMonitor(host_log, str(job_scope), sample_fn=enriched_snapshot)
            host.observe()
            host_observer = PeriodicHostObserver(host, host_interval_s)
            memory_before = periodic.read_job_memory(job_scope)
            save(directory / "job_memory_before.json", memory_before)
            points = DiskObservations(directory)
            if release_guard is not None:
                release_guard(directory)
                checked_at = time.monotonic()
            points.append(periodic.read_point(ready["pid"], ready["cgroup"], job_id,
                                             directory / "failed_point.json"))
            host_observer.start()
            if release_guard is not None and time.monotonic() - checked_at > 1:
                save(directory / "release_freshness_failed.json", dict(status="release_guard_stale"))
                raise ValueError("Release guard stale after initial observation")
            save(directory / "go.json", {"go": True})
            start, index = time.monotonic(), 0
            while True:
                index += 1
                time.sleep(max(0, start + index * interval_s - time.monotonic()))
                completed = (directory / "done.json").exists()
                if process.poll() is not None:
                    raise RuntimeError("Worker exited before final observation")
                if completed:
                    break
                if time.monotonic() - start > periodic.TIMEOUT + 30:
                    raise TimeoutError("Native command exceeded timeout and cleanup allowance")
            points.append(periodic.read_point(ready["pid"], ready["cgroup"], job_id,
                                             directory / "failed_point.json"))
            done = json.loads((directory / "done.json").read_text())
            status = "command_exited_zero" if done["exit_code"] == 0 else "command_failed"
            wall = (done["finished_ns"] - done["started_ns"]) / 1e9
            native_outcome(dict(status=status, native_wall_s=wall), done)
            completion = completion_evidence(points[0]["thread_affinity"], points[-1]["thread_affinity"],
                                             ready["pid"], done["finished_ns"])
            save(directory / "native_completion.json", completion)
            if completion["errors"]:
                raise ValueError("Native subtree completion unverified; retain failed attempt")
            host_summary = host_observer.finish(done["started_ns"] / 1e9, done["finished_ns"] / 1e9)
            save(directory / "host_process_summary.json", host_summary)
            memory = periodic.step_memory(periodic.interval_point(points[-1], job_id))
            save(directory / "step_memory.json", memory)
            memory_after = periodic.read_job_memory(job_scope)
            save(directory / "job_memory_after.json", memory_after)
            save(directory / "release.json", {"release": True})
            if process.wait(timeout=45) != 0:
                raise RuntimeError("Native step wrapper failed")
            with observe_reporting(directory, job_id, job_scope, periodic.read_job_memory):
                measured = dict(status=status, schema=SCHEMA, collector_arm="boundary",
                    policy=dict(POLICY), native=done, native_wall_s=wall, job_id=job_id,
                    launched=launched, placement=ready["placement"], point_records=points.records(),
                    native_completion=completion, step_memory=memory,
                    job_memory=dict(before=memory_before, after=memory_after),
                    host_process_observation=host_summary,
                    screening=periodic.evaluate_lineage(points, done, job_id),
                    root_context=periodic.evaluate(points, job_id),
                    scientific_timings_admitted=False, controlled_workload_verified=False,
                    native_outputs_validated=False, publication_ready=False,
                    limitations=["Two native boundary observations, not periodic native containment or affinity evidence.",
                        "Uses the unchanged periodic collector's worker, native clock, timeout and gate.",
                        "Common 30-second host observer is retained; its cost is not isolated by the paired contrast.",
                        "Screening flags and non-atomic windows are retained, not timing eligibility or quiet-host certification.",
                        "Native-step and job/reporting peaks overlap; no addition or baseline subtraction.",
                        "No pair-output equality, frozen-runtime/scheduler audit or engineering budget admission here."])
                save(directory / "boundary_report.json", measured)
            return measured
        except BaseException as error:
            save(directory / "boundary_failure.json", dict(status="boundary_collection_failed",
                error_type=type(error).__name__, error=str(error), scientific_timings_admitted=False))
            raise
        finally:
            try:
                if host_observer is not None:
                    host_observer.close()
            finally:
                if not (directory / "go.json").exists():
                    save(directory / "go.json", {"abort": True})
                if not (directory / "release.json").exists():
                    save(directory / "release.json", {"release": True})
                if process.poll() is None:
                    process.wait(timeout=periodic.TIMEOUT + 90)

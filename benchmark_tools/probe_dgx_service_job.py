"""Exercise Samwise suppression across one held/released five-second Slurm job."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import signal
import shlex
import sys
import time

from benchmark_tools.dgx_service_guard import ServiceGuard
from benchmark_tools.probe_dgx_step_separation import save


def wait_restore(guard, seconds=90):
    deadline = time.monotonic() + seconds
    while True:
        guard.check()
        try:
            guard.restore()
            return
        except ValueError as error:
            if "not authoritatively terminal" not in str(error):
                raise
        if time.monotonic() >= deadline:
            raise TimeoutError("Job not verified terminal; suppression retained")
        time.sleep(2)


def run(output, scaling_collector=False):
    guard = ServiceGuard(output)
    result = {"status": "prelaunch", "benchmark_submitted": False,
              "scientific_timings_admitted": False}
    paths = [Path(__file__).resolve().with_name(name) for name in (
        "probe_dgx_service_job.py", "dgx_service_guard.py",
        "probe_dgx_step_separation.py", "probe_host_counters.py")]
    if scaling_collector:
        paths = sorted(Path(__file__).resolve().parent.glob("*.py"))
    result["sources"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

    def interrupt(signum, frame):
        raise InterruptedError(f"Job probe signal {signum}")

    previous = {sig: signal.signal(sig, interrupt) for sig in (signal.SIGHUP, signal.SIGINT, signal.SIGTERM)}
    try:
        guard.begin()
        if guard.command(["squeue", "-h", "-p", "spark", "-o", "%i %T"]).strip():
            raise ValueError("Require empty spark queue before diagnostic submission")
        guard.before_submission()
        wrapped = "/bin/sleep 5"
        if scaling_collector:
            wrapped = shlex.join([sys.executable, "-B", "-m",
                "benchmark_tools.probe_scaling_collector", "--output", str(output / "collector")])
        raw = guard.command(["sbatch", "--hold", "--parsable", "--partition=spark",
            "--nodelist=spark-7ff0", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
            "--mem=96G", "--exclusive", "--time=" + ("00:03:00" if scaling_collector else "00:01:00"), "--no-requeue",
            "--job-name=service_guard_probe", "--output=" + str(output / "job.log"),
            "--chdir=" + str(Path(__file__).resolve().parent.parent), "--wrap=" + wrapped])
        match = re.fullmatch(r"([1-9][0-9]*)(?:;[A-Za-z0-9_.-]+)?\s*", raw)
        if match is None:
            raise ValueError("Submission identity unresolved; inspect scheduler without resubmitting")
        job = int(match[1])
        guard.bind_job(job)
        result["job_id"] = job
        try:
            guard.restore()
        except ValueError as error:
            if "not authoritatively terminal" not in str(error):
                raise
            guard.check()
            result["nonterminal_restoration_rejected"] = True
        else:
            raise RuntimeError("Held-job restoration unexpectedly succeeded")
        guard.command(["scontrol", "release", str(job)])
        wait_restore(guard, seconds=240 if scaling_collector else 90)
        if scaling_collector:
            receipt = json.loads((output / "collector" / "probe.json").read_text())
            if receipt["status"] != "collector_probe_completed" or receipt["job_id"] != job:
                raise ValueError("Collector probe did not complete for this job")
            result["collector"] = receipt
        result.update(status="job_lifecycle_probe_completed", restoration_verified=True)
    except BaseException as error:
        result.update(status="probe_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        try:
            if guard.stopped and not guard.restored:
                if guard.job is not None:
                    guard.command(["scancel", str(guard.job)])
                    wait_restore(guard)
                else:
                    guard.restore()
            result["restoration_verified"] = guard.restored
            if result["sources"] != {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}:
                raise ValueError("Job probe source changed")
        except BaseException as error:
            result.update(cleanup_error_type=type(error).__name__, cleanup_error=str(error))
            raise
        finally:
            for sig, handler in previous.items():
                signal.signal(sig, handler)
            if output.is_dir():
                save(output / "probe.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scaling-collector", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run(args.output.resolve(), args.scaling_collector), sort_keys=True))

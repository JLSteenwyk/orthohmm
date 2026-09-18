"""Counter-only native-step collector for engineering smokes, not timing admission."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe_host_counters import snapshot, summarize
from probe_dgx_step_separation import save, wait_file, validate_scopes


def validate(command, cpus, timeout_s, interval_s):
    if not command or not all(isinstance(x, str) and x for x in command) or not Path(command[0]).is_absolute():
        raise ValueError("Require absolute executable and nonempty arguments")
    if cpus != 20 or not 0 < timeout_s <= 900 or interval_s != 1.:
        raise ValueError("Require frozen 20CPU/900s-maximum/1s smoke limits")


def worker(directory):
    plan = json.loads((directory / "command.json").read_text())
    validate(plan["command"], plan["cpus"], plan["timeout_s"], plan["interval_s"])
    save(directory / "native_before.json", snapshot())
    save(directory / "ready.json", {"pid": os.getpid()})
    wait_file(directory / "go.json")
    started = time.monotonic_ns()
    with (directory / "native.log").open("x") as log:
        try:
            completed = subprocess.run(plan["command"], stdout=log, stderr=subprocess.STDOUT,
                                       timeout=plan["timeout_s"])
            code = completed.returncode
        except subprocess.TimeoutExpired:
            # The owning Slurm allocation bounds descendants if the direct child times out.
            code = 124
    finished = time.monotonic_ns()
    save(directory / "native_after.json", snapshot())
    save(directory / "done.json", dict(exit_code=code, started_ns=started, finished_ns=finished))


def measure(command, directory, job_id, cpus, memory_bytes, timeout_s, interval_s,
            monitor_host=True, host_interval_s=30.):
    validate(command, cpus, timeout_s, interval_s)
    if int(os.environ["SLURM_JOB_ID"]) != job_id or os.environ.get("SLURM_CPUS_PER_TASK") != "20":
        raise ValueError("Wrong scheduler allocation")
    if os.uname().nodename != "spark-7ff0" or memory_bytes != 96 * 1024 ** 3:
        raise ValueError("Wrong host or memory request")
    directory = Path(directory).absolute()
    directory.mkdir(exist_ok=False)
    save(directory / "command.json", dict(command=command, cpus=cpus, timeout_s=timeout_s, interval_s=interval_s))
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
                sys.executable, "-B", str(Path(__file__).resolve()), "--worker", str(directory)]
    with (directory / "step.log").open("x") as log:
        process = subprocess.Popen(launched, stdout=log, stderr=subprocess.STDOUT)
        try:
            wait_file(directory / "ready.json")
            before = json.loads((directory / "native_before.json").read_text())
            samples = [snapshot()]
            validate_scopes(samples[0], before, job_id)
            save(directory / "go.json", {"go": True})
            deadline = time.monotonic() + timeout_s + 60
            with (directory / "observer.jsonl").open("x") as stream:
                stream.write(json.dumps(samples[0]) + "\n")
                while process.poll() is None:
                    if time.monotonic() > deadline:
                        raise TimeoutError("Native step exceeded command timeout and cleanup allowance")
                    time.sleep(interval_s)
                    sample = snapshot()
                    samples.append(sample)
                    stream.write(json.dumps(sample) + "\n")
                    stream.flush()
            if process.returncode != 0:
                raise RuntimeError("Native step wrapper failed")
            done = json.loads((directory / "done.json").read_text())
            after = json.loads((directory / "native_after.json").read_text())
            for sample in samples:
                validate_scopes(sample, before, job_id)
            summary = summarize(before, after, os.sysconf("SC_CLK_TCK"))
            if not before["finished_monotonic_ns"] < done["started_ns"] < done["finished_ns"] < after["started_monotonic_ns"]:
                raise ValueError("Native work not bracketed")
            errors = sum(len(s["errors"]) for s in [before, after, *samples])
            result = dict(status="command_exited_zero" if done["exit_code"] == 0 else "command_failed",
                native=done, native_wall_s=(done["finished_ns"] - done["started_ns"]) / 1e9,
                job_id=job_id, launched=launched, native_snapshots=[before, after], host_summary=summary,
                observer_snapshots=len(samples), counter_read_errors=errors, publication_ready=False,
                scientific_results_admitted=False, controlled_workload_verified=False,
                limitations=["Engineering smoke; no exclusivity or calibrated overhead claim.",
                    "Observer shares allocation CPUs with the native step; cgroup scopes differ.",
                    "Cgroup memory is not RSS and includes the native wrapper.",
                    "Host and native counters have different read windows; no foreign-CPU subtraction.",
                    "Native outputs require separate validation; nonzero/read-error runs are retained."])
            save(directory / "counter_report.json", result)
            return result
        finally:
            if process.poll() is None:
                if not (directory / "go.json").exists():
                    save(directory / "go.json", {"go": True, "failed_observation": True})
                process.wait(timeout=timeout_s + 90)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, required=True)
    worker(parser.parse_args().worker)

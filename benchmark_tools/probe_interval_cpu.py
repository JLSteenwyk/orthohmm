"""Fixed-window interval CPU controls; not scientific timing admission."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.probe_host_counters import snapshot, summarize, parse_group
from benchmark_tools.probe_dgx_step_separation import save, wait_file, validate_scopes, scope_parts, burn
from benchmark_tools.screen_bracketed_cpu import usage

INTERVALS = 10
PERIOD_S = 1.
MAX_GAP_S = 1.5


def cpu_usage(point):
    return usage({"optional": {"cgroup_cpu.stat": point["native_cpu_stat"]}})


def validate_point(point, job):
    before, after = point["host"]
    if before["errors"] or after["errors"]:
        raise ValueError("Counter read failures")
    summarize(before, after, point["ticks_per_second"])
    native = {"raw": {"cgroup_membership": point["native_membership"]}}
    validate_scopes(before, native, job)
    validate_scopes(after, native, job)
    expected = str(Path(*scope_parts(native, job)))
    if point["native_cpu_scope"] != expected:
        raise ValueError("CPU counter is not the native step aggregate")
    start, end = point["native_read_ns"]
    if (type(start) is not int or type(end) is not int
            or not before["finished_monotonic_ns"] <= start <= end <= after["started_monotonic_ns"]):
        raise ValueError("Host reads do not enclose native read")
    cpu_usage(point)


def interval(left, right, job, enforce_gap=True):
    for point in (left, right):
        validate_point(point, job)
    for key in ("native_membership", "native_cpu_scope", "ticks_per_second"):
        if left[key] != right[key]:
            raise ValueError("Changed interval identity")
    if left["host"][1]["finished_monotonic_ns"] >= right["host"][0]["started_monotonic_ns"]:
        raise ValueError("Overlapping observation points")
    duration = (sum(right["native_read_ns"]) - sum(left["native_read_ns"])) / 2e9
    if enforce_gap and not .5 <= duration <= MAX_GAP_S:
        raise ValueError("Missing or irregular observation interval")
    host = summarize(left["host"][0], right["host"][1], left["ticks_per_second"])
    usec = cpu_usage(right) - cpu_usage(left)
    if usec < 0:
        raise ValueError("Native CPU counter decreased")
    residual = host["accounted_host_busy_cpu_s"] - usec / 1e6
    reasons = []
    if residual / duration > .25:
        reasons.append("excess_unassigned_cpu")
    if residual < -.5:
        reasons.append("negative_accounting_discrepancy")
    if host["cpu_delta_ticks"]["steal"]:
        reasons.append("host_steal_time")
    return dict(wall_s=duration, native_cpu_s=usec / 1e6,
        host_busy_cpu_s=host["accounted_host_busy_cpu_s"], signed_unassigned_cpu_s=residual,
        signed_unassigned_average_cores=residual / duration, reasons=reasons,
        screen_passed=not reasons,
        outer_read_overhang_s=((left["native_read_ns"][0] - left["host"][0]["started_monotonic_ns"])
            + (right["host"][1]["finished_monotonic_ns"] - right["native_read_ns"][1])) / 1e9)


def evaluate(points, job):
    if len(points) != INTERVALS + 1:
        raise ValueError("Incomplete fixed-window observation")
    intervals = [interval(a, b, job) for a, b in zip(points, points[1:])]
    whole = interval(points[0], points[-1], job, enforce_gap=False)
    return dict(intervals=intervals, whole_window=whole,
        flagged_intervals=[i for i, value in enumerate(intervals) if not value["screen_passed"]],
        controlled_workload_verified=False, scientific_timings_admitted=False)


def read_point(ready, job):
    pid = ready["pid"]
    membership = Path(f"/proc/{pid}/cgroup").read_text()
    if membership != ready["cgroup"]:
        raise ValueError("Worker scope changed")
    native = {"raw": {"cgroup_membership": membership}}
    scope = str(Path(*scope_parts(native, job)))
    before = snapshot()
    start = time.monotonic_ns()
    raw = (Path("/sys/fs/cgroup") / scope.lstrip("/") / "cpu.stat").read_text()
    end = time.monotonic_ns()
    after = snapshot()
    if Path(f"/proc/{pid}/cgroup").read_text() != membership:
        raise ValueError("Worker disappeared or changed scope")
    point = dict(host=[before, after], native_membership=membership, native_cpu_scope=scope,
        native_read_ns=[start, end], native_cpu_stat=raw, ticks_per_second=os.sysconf("SC_CLK_TCK"))
    validate_point(point, job)
    return point


def worker(directory):
    save(directory / "ready.json", dict(pid=os.getpid(), cgroup=Path("/proc/self/cgroup").read_text()))
    wait_file(directory / "go.json")
    started, cpu_start = time.monotonic_ns(), time.process_time()
    value = 1
    while not (directory / "release.json").exists():
        if (time.monotonic_ns() - started) / 1e9 > 40:
            raise TimeoutError("Worker release missing")
        for _ in range(10000):
            value = (value * 1664525 + 1013904223) % 4294967296
    save(directory / "done.json", dict(started_ns=started, finished_ns=time.monotonic_ns(),
        cpu_s=time.process_time() - cpu_start, checksum=value))


def run(output):
    job = int(os.environ["SLURM_JOB_ID"])
    if os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "2":
        raise ValueError("Require DGX two-CPU engineering allocation")
    output.mkdir(exist_ok=False)
    source = Path(__file__).resolve()
    names = (source.name, "probe_host_counters.py", "probe_dgx_step_separation.py", "screen_bracketed_cpu.py")
    hashes = {name: hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in names}
    trials = []
    for mode in ("quiet", "completed_burst"):
        directory = output / mode
        directory.mkdir()
        command = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=1",
                   sys.executable, "-I", "-B", str(source), "--worker", str(directory)]
        with (directory / "step.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            sibling = None
            try:
                ready = wait_file(directory / "ready.json")
                points = [read_point(ready, job)]
                save(directory / "point_00.json", points[0])
                save(directory / "go.json", {"go": True})
                start = time.monotonic()
                for index in range(1, INTERVALS + 1):
                    time.sleep(max(0, start + index * PERIOD_S - time.monotonic()))
                    points.append(read_point(ready, job))
                    save(directory / f"point_{index:02d}.json", points[-1])
                    if index == 2 and mode == "completed_burst":
                        sibling = subprocess.Popen([sys.executable, "-I", "-B", str(source), "--burn"],
                                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                save(directory / "release.json", {"release": True})
                if process.wait(timeout=45) != 0:
                    raise RuntimeError("Native control failed")
                load = None
                if sibling is not None:
                    stdout, stderr = sibling.communicate(timeout=10)
                    if sibling.returncode:
                        raise RuntimeError(stderr)
                    load = json.loads(stdout)
                    if not .75 <= load["cpu_seconds"] < .85:
                        raise ValueError("Unexpected burst CPU duration")
                    if parse_group(load["cgroup"]) != parse_group(points[0]["host"][0]["raw"]["cgroup_membership"]):
                        raise ValueError("Burst escaped observer scope")
                    if not (points[2]["host"][1]["finished_monotonic_ns"] < load["started_ns"]
                            < load["finished_ns"] < points[4]["host"][0]["started_monotonic_ns"]):
                        raise ValueError("Burst did not finish within frozen observation window")
                result = evaluate(points, job)
                met = (not result["flagged_intervals"] if mode == "quiet" else bool(result["flagged_intervals"]))
                met = met and result["whole_window"]["screen_passed"]
                trial = dict(mode=mode, ready=ready, command=command, points=points, sibling=load,
                    native=json.loads((directory / "done.json").read_text()), result=result,
                    expectation_met=met)
                save(directory / "trial.json", trial)
                trials.append(trial)
            finally:
                if not (directory / "release.json").exists():
                    save(directory / "release.json", {"release": True, "failed_observation": True})
                if not (directory / "go.json").exists():
                    save(directory / "go.json", {"go": True, "failed_observation": True})
                process.wait(timeout=45)
                if sibling is not None and sibling.poll() is None:
                    sibling.communicate(timeout=15)
    if hashes != {name: hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in names}:
        raise ValueError("Probe sources changed")
    save(output / "report.json", dict(status="interval_cpu_controls_evaluated", job_id=job,
        sources=hashes, trials=trials, all_control_expectations_met=all(t["expectation_met"] for t in trials),
        controlled_workload_verified=False, scientific_timings_admitted=False, publication_ready=False,
        limitations=["Two fixed-order engineering controls, not a scientific timing panel.",
            "Ten-second observation window excludes some worker startup and teardown.",
            "Intervals overlap in outer host read brackets; residuals must not be summed.",
            "Residual includes observer and kernel work; no identified foreign-load or rigorous bound.",
            "Accounting delay, non-CPU interference and native-pipeline overhead remain unresolved."]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--worker", type=Path)
    mode.add_argument("--burn", action="store_true")
    mode.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker)
    elif args.burn:
        started = time.monotonic_ns()
        result = burn()
        result.update(started_ns=started, finished_ns=time.monotonic_ns())
        print(json.dumps(result))
    else:
        run(args.output.absolute())

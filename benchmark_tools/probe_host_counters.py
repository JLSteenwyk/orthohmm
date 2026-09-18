"""Read-only cumulative host counters for prospective observer development."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import platform
import time

CPU_FIELDS = ("user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal", "guest", "guest_nice")
BUSY_FIELDS = ("user", "nice", "system", "irq", "softirq")


def parse_cpu(text):
    rows = [line.split() for line in text.splitlines() if line.startswith("cpu ")]
    if len(rows) != 1 or len(rows[0]) != 11:
        raise ValueError("Require one ten-field aggregate CPU row")
    result = dict(zip(CPU_FIELDS, map(int, rows[0][1:])))
    if any(v < 0 for v in result.values()):
        raise ValueError("Negative CPU counter")
    return result


def parse_group(text):
    rows = [line[3:] for line in text.splitlines() if line.startswith("0::")]
    if len(rows) != 1:
        raise ValueError("Require unified cgroup membership")
    path = PurePosixPath(rows[0])
    if not path.is_absolute() or ".." in path.parts:
        raise ValueError("Invalid cgroup path")
    return str(path)


def snapshot():
    started = time.monotonic_ns()
    membership = Path("/proc/self/cgroup").read_text()
    scope = parse_group(membership)
    group = Path("/sys/fs/cgroup") / scope.lstrip("/")
    raw = {"proc_stat": Path("/proc/stat").read_text(),
           "online_cpus": Path("/sys/devices/system/cpu/online").read_text(),
           "boot_id": Path("/proc/sys/kernel/random/boot_id").read_text(),
           "cgroup_membership": membership}
    optional, errors = {}, []
    paths = {**{f"host_{resource}_pressure": Path("/proc/pressure") / resource for resource in ("cpu", "memory", "io")},
             **{f"cgroup_{name}": group / name for name in ("cpu.stat", "memory.current", "memory.peak", "memory.events")}}
    for name, path in paths.items():
        try:
            optional[name] = path.read_text()
        except OSError as error:
            errors.append({"field": name, "type": type(error).__name__, "errno": error.errno})
    if parse_group(Path("/proc/self/cgroup").read_text()) != scope:
        raise ValueError("Observer cgroup changed")
    return dict(started_monotonic_ns=started, finished_monotonic_ns=time.monotonic_ns(),
                cpu_ticks=parse_cpu(raw["proc_stat"]), raw=raw, optional=optional, errors=errors)


def summarize(before, after, ticks_per_second):
    if type(ticks_per_second) is not int or ticks_per_second <= 0:
        raise ValueError("Invalid clock tick rate")
    for key in ("boot_id", "online_cpus", "cgroup_membership"):
        if before["raw"][key] != after["raw"][key]:
            raise ValueError("Changed host or observer scope: " + key)
    times = [before["started_monotonic_ns"], before["finished_monotonic_ns"],
             after["started_monotonic_ns"], after["finished_monotonic_ns"]]
    if any(type(t) is not int or t < 0 for t in times) or not times[0] <= times[1] < times[2] <= times[3]:
        raise ValueError("Invalid or overlapping snapshot times")
    for sample in (before, after):
        if sample["cpu_ticks"] != parse_cpu(sample["raw"]["proc_stat"]):
            raise ValueError("Parsed counters disagree with raw evidence")
    delta = {k: after["cpu_ticks"][k] - before["cpu_ticks"][k] for k in CPU_FIELDS}
    if any(v < 0 for k, v in delta.items() if k != "iowait"):
        raise ValueError("CPU counter decreased")
    span = ((times[2] + times[3]) - (times[0] + times[1])) / 2e9
    busy = sum(delta[k] for k in BUSY_FIELDS) / ticks_per_second
    return dict(status="host_counter_diagnostic_only", cpu_delta_ticks=delta,
        midpoint_span_s=span, accounted_host_busy_cpu_s=busy, accounted_host_busy_average_cores=busy / span,
        iowait_decreased=delta["iowait"] < 0, controlled_workload_verified=False)


def run(output, interval):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if not math.isfinite(interval) or not 0 < interval <= 60:
        raise ValueError("Require interval in (0, 60] seconds")
    source = Path(__file__).resolve()
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    before = snapshot()
    time.sleep(interval)
    after = snapshot()
    ticks = os.sysconf("SC_CLK_TCK")
    result = dict(status="prospective_host_counter_probe", hostname=platform.node(), kernel=platform.release(),
        source=dict(path=str(source), sha256=digest), python_version=platform.python_version(),
        clock_ticks_per_second=ticks, requested_interval_s=interval, snapshots=[before, after],
        summary=summarize(before, after, ticks), publication_ready=False,
        limitations=["Read-only availability probe, not a timing run, exclusivity proof or observer calibration.",
            "Host counters include this observer, SSH and unrelated workloads; no foreign-CPU subtraction is made.",
            "Guest time is already included in user/nice and is not added again; steal is reported separately.",
            "Reads are non-atomic; midpoint rates are descriptive, not rigorous contention bounds.",
            "Optional fields retain read failures. Pressure counters do not prove absence of interference.",
            "Cannot retroactively supply missing evidence for previous timing runs."])
    if hashlib.sha256(source.read_bytes()).hexdigest() != digest:
        raise ValueError("Probe source changed")
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--interval", type=float, default=3.)
    args = parser.parse_args()
    run(args.output.absolute(), args.interval)

"""Bounded read-only DGX service/device evidence, never an isolation certificate."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


COMMANDS = {
    "system_services": ["systemctl", "list-units", "--type=service", "--state=running", "--no-pager", "--plain", "--no-legend"],
    "user_services": ["systemctl", "--user", "list-units", "--type=service", "--state=running", "--no-pager", "--plain", "--no-legend"],
    "timers": ["systemctl", "list-timers", "--all", "--no-pager", "--no-legend"],
    "containers": ["docker", "ps", "--format", "{{.ID}} {{.Image}} {{.Status}}"],
    "ollama_models": ["ollama", "ps"],
    "gpu": ["nvidia-smi", "--query-gpu=name,uuid,driver_version,utilization.gpu,utilization.memory,memory.used,temperature.gpu,power.draw", "--format=csv,noheader,nounits"],
    "gpu_compute": ["nvidia-smi", "--query-compute-apps=pid,process_name,used_memory", "--format=csv,noheader,nounits"],
    "scheduler": ["squeue", "-h", "-p", "spark", "-o", "%i %T %N"],
}
STATIC_COMMANDS = {"cpu": ["lscpu", "-J"], "kernel": ["uname", "-a"]}


def collect(static=False):
    if os.uname().nodename != "spark-7ff0":
        raise ValueError("Require the designated DGX host")
    result = dict(started_unix_ns=time.time_ns(), started_monotonic_ns=time.monotonic_ns(),
        host=os.uname().nodename, commands={}, files={}, static=static,
        environmental_validity_established=False,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    for key, argv in (COMMANDS | (STATIC_COMMANDS if static else {})).items():
        row = dict(command=argv, started_monotonic_ns=time.monotonic_ns())
        try:
            completed = subprocess.run(argv, capture_output=True, text=True, timeout=10,
                                       env=dict(os.environ, LC_ALL="C", SYSTEMD_COLORS="0"))
            row.update(returncode=completed.returncode, stdout=completed.stdout, stderr=completed.stderr)
        except (OSError, subprocess.TimeoutExpired) as error:
            row.update(error_type=type(error).__name__, error=str(error))
        row["finished_monotonic_ns"] = time.monotonic_ns()
        result["commands"][key] = row
    paths = [Path(p) for p in ("/proc/diskstats", "/proc/meminfo", "/proc/pressure/cpu",
        "/proc/pressure/io", "/proc/pressure/memory", "/proc/sys/kernel/random/boot_id")]
    paths.extend(Path("/sys/class/thermal").glob("thermal_zone*/temp"))
    paths.extend(Path("/sys/devices/system/cpu/cpufreq").glob("policy*/scaling_cur_freq"))
    for path in paths:
        row = dict(started_monotonic_ns=time.monotonic_ns())
        try:
            row["text"] = path.read_text()
        except OSError as error:
            row.update(error_type=type(error).__name__, error=str(error))
        row["finished_monotonic_ns"] = time.monotonic_ns()
        result["files"][str(path)] = row
    result.update(finished_unix_ns=time.time_ns(), finished_monotonic_ns=time.monotonic_ns(),
        limitations=["Sequential snapshots are not atomic and miss activity between observations.",
            "Device counters do not attribute I/O to the measured job; GPU N/A is not zero.",
            "Service presence is not permission for arbitrary work; configuration policy is separate."])
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--static", action="store_true")
    args = parser.parse_args()
    with args.output.open("x") as handle:
        json.dump(collect(args.static), handle, indent=2, sort_keys=True)
        handle.write("\n")

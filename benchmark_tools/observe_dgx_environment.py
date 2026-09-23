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
    "system_configuration": ["systemctl", "show", "*.service", "*.timer",
        "--property=Id,LoadState,FragmentPath,DropInPaths,NeedDaemonReload"],
    "user_configuration": ["systemctl", "--user", "show", "*.service", "*.timer",
        "--property=Id,LoadState,FragmentPath,DropInPaths,NeedDaemonReload"],
}
STATIC_COMMANDS = {"cpu": ["lscpu", "-J"], "kernel": ["uname", "-a"]}


def configuration_fingerprints(raw):
    units = {}
    for block in raw.strip().split("\n\n"):
        fields = {}
        for line in block.splitlines():
            key, separator, value = line.partition("=")
            if not separator or key in fields:
                raise ValueError("Malformed or duplicate unit property")
            fields[key] = value
        if set(fields) != {"Id", "LoadState", "FragmentPath", "DropInPaths", "NeedDaemonReload"}:
            raise ValueError("Incomplete unit configuration properties")
        name = fields["Id"]
        if not name.endswith((".service", ".timer")) or name in units:
            raise ValueError("Invalid or duplicate unit identity")
        files = []
        for value in [fields["FragmentPath"], *fields["DropInPaths"].split()]:
            if not value:
                continue
            row = {"path": value}
            try:
                path = Path(value)
                if not path.is_absolute() or any(c in value for c in ('\\', '"', "'")):
                    raise ValueError("Unsupported escaped unit path")
                row["resolved_path"] = str(path.resolve(strict=True))
                row["symlink_target"] = os.readlink(path) if path.is_symlink() else None
                if not path.is_file() or path.stat().st_size > 1024**2:
                    raise ValueError("Unit path is not a bounded regular file")
                data = path.read_bytes()
                row.update(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
            except (OSError, ValueError) as error:
                row.update(error_type=type(error).__name__, error=str(error))
            files.append(row)
        units[name] = {"properties": fields, "files": files}
    return units


class Recorder:
    """Stream session evidence without blocking the native one-second collector."""

    def __init__(self, directory):
        self.directory = Path(directory)
        self.index = 0
        self.next_due = 0.

    def observe(self, stage, static=False):
        self.directory.mkdir(exist_ok=True)
        result = collect(static)
        result.update(session_stage=stage, index=self.index)
        with (self.directory / f"snapshot_{self.index:06d}.json").open("x") as handle:
            json.dump(result, handle, indent=2, sort_keys=True)
            handle.write("\n")
        self.index += 1
        self.next_due = time.monotonic() + 30.

    def periodic(self):
        if time.monotonic() >= self.next_due:
            self.observe("job_wait")


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
    result["unit_configurations"] = {}
    for key in ("system_configuration", "user_configuration"):
        command = result["commands"][key]
        try:
            if command.get("returncode") != 0:
                raise ValueError("Configuration command did not succeed")
            result["unit_configurations"][key] = configuration_fingerprints(command["stdout"])
        except ValueError as error:
            result["unit_configurations"][key] = {"error_type": type(error).__name__, "error": str(error)}
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
            "Service presence is not permission for arbitrary work; configuration policy is separate.",
            "Unit fingerprints cover loaded service/timer fragments and drop-ins, not all external service configuration."])
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--static", action="store_true")
    args = parser.parse_args()
    with args.output.open("x") as handle:
        json.dump(collect(args.static), handle, indent=2, sort_keys=True)
        handle.write("\n")

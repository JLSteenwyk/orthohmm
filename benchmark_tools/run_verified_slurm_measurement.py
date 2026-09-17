"""Bracket the unchanged Slurm collector with separately timed runtime checks."""

import argparse
import importlib
import json
import os
from pathlib import Path
import sys
import time

from snapshot_runtime_trees import digest, verify


def check_manifests(specifications):
    checks = []
    for path, expected_sha in specifications:
        path = Path(path)
        if digest(path) != expected_sha:
            raise ValueError("Runtime manifest digest differs: " + str(path))
        expected = json.loads(path.read_text())
        checked = verify(expected)
        if digest(path) != expected_sha:
            raise ValueError("Runtime manifest changed while checking")
        checks.append({"path": str(path), "sha256": expected_sha, **checked})
    if not checks:
        raise ValueError("Require at least one pinned runtime manifest")
    return checks


def run_checked(specifications, output, measurement, checker=check_manifests):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    result = {"status": "preflight", "scientific_results_admitted": False,
              "source_sha256": digest(Path(__file__)),
              "limitations": ["Identity checks are outside the native command timer.",
                              "Hashing warms file caches; this is not cold-cache timing.",
                              "Cgroup memory peak may retain preparation and verification allocations.",
                              "Native outputs and collector records still require independent admission."]}
    before_ok = False
    try:
        start = time.monotonic()
        result["before"] = checker(specifications)
        result["before_check_wall_s"] = time.monotonic() - start
        before_ok = True
        result["measurement"] = measurement(output / "measurement")
        result["status"] = result["measurement"]["status"]
    except Exception as error:
        result.update(status="verified_wrapper_failed", error_type=type(error).__name__, error=str(error))
    finally:
        if before_ok:
            start = time.monotonic()
            try:
                result["after"] = checker(specifications)
            except Exception as error:
                result.update(status="runtime_changed_or_unverifiable", after_error_type=type(error).__name__,
                              after_error=str(error))
            result["after_check_wall_s"] = time.monotonic() - start
        with (output / "verification.json").open("x") as handle:
            json.dump(result, handle, indent=2, sort_keys=True)
            handle.write("\n")
    return result


def load_collector(directory, specifications):
    source = directory / "measure_slurm_command.py"
    rows = [row for path, _ in specifications for row in json.loads(Path(path).read_text())["records"]]
    if not any(row["path"] == str(source) and row.get("sha256") == digest(source) for row in rows):
        raise ValueError("Collector source is not in the verified runtime inventory")
    sys.path.insert(0, str(directory))
    module = importlib.import_module("measure_slurm_command")
    if Path(module.__file__).resolve() != source.resolve():
        raise ValueError("Wrong collector import path")
    return module.measure


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime", nargs=2, action="append", required=True, metavar=("MANIFEST", "SHA256"))
    for flag in ("collector", "output", "cwd"):
        parser.add_argument("--" + flag, type=Path, required=True)
    for flag in ("job-id", "cpus", "memory-gib"):
        parser.add_argument("--" + flag, type=int, required=True)
    parser.add_argument("--timeout", type=float, required=True)
    parser.add_argument("--interval", type=float, default=1.)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or not Path(command[0]).is_absolute():
        raise ValueError("Require explicit absolute native executable")
    if not sys.dont_write_bytecode or not sys.pycache_prefix or Path(sys.pycache_prefix).exists():
        raise ValueError("Require disabled bytecode writes and absent isolated cache prefix at startup")
    for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        if os.environ.get(key):
            raise ValueError("Unexpected dynamic loader override: " + key)
    if Path("/etc/ld.so.preload").exists():
        raise ValueError("System preload configuration requires separate review")
    os.chdir(args.cwd.resolve())

    def measurement(output):
        measure = load_collector(args.collector.resolve(), args.runtime)
        return measure(command, output, args.job_id, args.cpus, args.memory_gib * 1024 ** 3,
                       args.timeout, args.interval, monitor_host=True, host_interval_s=30.)

    result = run_checked(args.runtime, args.output.resolve(), measurement)
    raise SystemExit(0 if result["status"] == "command_exited_zero" else 1)

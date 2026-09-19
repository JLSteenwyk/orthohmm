"""Retain controller records before expiry; never admit timing measurements."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import time


TERMINAL = {"COMPLETED", "FAILED", "TIMEOUT", "OUT_OF_MEMORY", "CANCELLED",
            "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED"}
REQUIRED = {"JobId", "ArrayJobId", "ArrayTaskId", "JobState", "ExitCode",
            "Restarts", "Requeue", "NodeList", "OverSubscribe", "MinMemoryNode",
            "NumNodes", "NumCPUs", "CPUs/Task", "Command"}


def terminal_records(raw, array_id, tasks):
    """Select exact task records, never compressed pending-array placeholders."""
    records = {}
    for line in raw.splitlines():
        fields = {}
        for key, value in re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", line):
            if key in fields:
                raise ValueError("Duplicate scheduler field: " + key)
            fields[key] = value
        if fields.get("ArrayJobId") != str(array_id):
            continue
        task = fields.get("ArrayTaskId", "")
        if not task.isdigit():
            continue
        index = int(task)
        if index not in tasks or fields.get("JobState") not in TERMINAL:
            continue
        if (not REQUIRED <= fields.keys() or not fields["JobId"].isdigit()
                or int(fields["JobId"]) <= 0
                or not re.fullmatch(r"\d+:\d+", fields["ExitCode"])):
            raise ValueError("Incomplete detailed terminal record")
        if index in records:
            raise ValueError("Duplicate terminal task")
        records[index] = line + "\n"
    return records


def capture(array_id, count, output, *, interval=5, max_seconds=86400,
            run=subprocess.run, clock=time.monotonic, sleep=time.sleep):
    if type(array_id) is not int or array_id <= 0 or type(count) is not int or count <= 0:
        raise ValueError("Require positive integer array ID and task count")
    if any(not math.isfinite(v) or v <= 0 for v in (interval, max_seconds)):
        raise ValueError("Require positive finite polling interval and deadline")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    tasks = set(range(count))
    retained = {}
    errors = 0
    started = clock()
    poll = 0
    # Each raw observation is immutable; a missing later controller record
    # cannot erase a previously captured terminal allocation.
    while clock() - started < max_seconds and set(retained) != tasks:
        argv = ["scontrol", "show", "job", str(array_id), "--oneliner"]
        observation = {"argv": argv, "observed_unix_ns": time.time_ns()}
        try:
            result = run(argv, capture_output=True, text=True, timeout=10, check=False)
            observation.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
            if result.returncode:
                errors += 1
            else:
                try:
                    records = terminal_records(result.stdout, array_id, tasks)
                    for index, raw in records.items():
                        if index in retained:
                            continue
                        path = output / f"scheduler_{index}.txt"
                        with path.open("x") as handle:
                            handle.write(raw)
                        retained[index] = {"path": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                                           "poll": poll}
                except ValueError as exc:
                    observation["parse_error"] = str(exc)
                    errors += 1
        except (OSError, subprocess.TimeoutExpired) as exc:
            observation.update(error_type=type(exc).__name__, error=str(exc))
            errors += 1
        with (output / f"poll_{poll:06d}.json").open("x") as handle:
            json.dump(observation, handle, indent=2, sort_keys=True)
            handle.write("\n")
        poll += 1
        if set(retained) != tasks:
            sleep(min(interval, max(0, max_seconds - (clock() - started))))
    report = dict(status="complete_controller_capture" if set(retained) == tasks else "incomplete_controller_capture",
                  array_id=array_id, expected_tasks=count, retained=retained,
                  missing_tasks=sorted(tasks - set(retained)), polls=poll, observation_errors=errors,
                  scientific_timings_admitted=False,
                  limitations=["Collection only, not allocation-policy or native-output validation.",
                               "No recovery of records expired before observation; no sacct substitution.",
                               "Run on the controller host; no SSH or compute-node reads are performed.",
                               "Retains the first detailed terminal record; retries/requeues require separate policy validation."])
    with (output / "capture.json").open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--array-id", required=True, type=int)
    parser.add_argument("--tasks", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--interval", type=float, default=5)
    parser.add_argument("--max-seconds", type=float, default=86400)
    args = parser.parse_args()
    result = capture(args.array_id, args.tasks, args.output, interval=args.interval, max_seconds=args.max_seconds)
    return 0 if not result["missing_tasks"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

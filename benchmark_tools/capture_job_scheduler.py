"""Retain terminal controller records for an explicit non-array job list."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

from benchmark_tools.capture_array_scheduler import TERMINAL, REQUIRED


def terminal_record(raw, job):
    for line in raw.splitlines():
        pairs = re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", line)
        fields = dict(pairs)
        if len(fields) != len(pairs):
            raise ValueError("Duplicate scheduler field")
        if fields.get("JobId") != str(job):
            continue
        if "ArrayJobId" in fields or "ArrayTaskId" in fields:
            raise ValueError("Require non-array job")
        if fields.get("JobState") not in TERMINAL:
            return None
        if not REQUIRED - {"ArrayJobId", "ArrayTaskId"} <= fields.keys():
            raise ValueError("Incomplete terminal allocation")
        if not re.fullmatch(r"\d+:\d+", fields["ExitCode"]):
            raise ValueError("Invalid exit status")
        return line + "\n"
    return None


def capture(jobs, output, max_seconds=14400):
    if not jobs or len(set(jobs)) != len(jobs) or any(type(j) is not int or j <= 0 for j in jobs):
        raise ValueError("Require unique positive job IDs")
    output.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve()
    sources = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
               (source, source.with_name("capture_array_scheduler.py"))}
    start, poll, retained, errors = time.monotonic(), 0, {}, 0
    while len(retained) < len(jobs) and time.monotonic() - start < max_seconds:
        for job in jobs:
            if job in retained:
                continue
            observation = dict(job_id=job, observed_unix_ns=time.time_ns())
            try:
                p = subprocess.run(["scontrol", "show", "job", str(job), "--oneliner"],
                                   capture_output=True, text=True, timeout=10)
                observation.update(returncode=p.returncode, stdout=p.stdout, stderr=p.stderr)
                if p.returncode:
                    errors += 1
                else:
                    raw = terminal_record(p.stdout, job)
                    if raw is not None:
                        path = output / f"scheduler_{job}.txt"
                        with path.open("x") as stream:
                            stream.write(raw)
                        retained[job] = dict(path=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            except (OSError, ValueError, subprocess.TimeoutExpired) as error:
                observation.update(error_type=type(error).__name__, error=str(error))
                errors += 1
            with (output / f"poll_{poll:06d}_{job}.json").open("x") as stream:
                json.dump(observation, stream, indent=2, sort_keys=True)
                stream.write("\n")
        poll += 1
        if len(retained) < len(jobs):
            time.sleep(5)
    for name, digest in sources.items():
        if hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() != digest:
            raise ValueError("Recorder source changed")
    result = dict(status="complete" if len(retained) == len(jobs) else "incomplete", jobs=jobs,
                  retained=retained, missing=[j for j in jobs if j not in retained], polls=poll,
                  observation_errors=errors, sources=sources, scientific_timings_admitted=False)
    with (output / "capture.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", nargs="+", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(0 if capture(args.jobs, args.output)["status"] == "complete" else 1)

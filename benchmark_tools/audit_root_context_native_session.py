"""Bind the one-hour native panel waiting receipts to terminal Slurm evidence."""

import argparse
from datetime import datetime
import json
from pathlib import Path
import re
from zoneinfo import ZoneInfo

from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.submit_root_context_native_session import (
    command, job_id, QUEUE_COMMAND, LOCAL_TIMEOUT, REMOTE_TIMEOUT, KILL_GRACE, RECIPE_ROOT)
from benchmark_tools.verify_root_context_native_provenance import scheduler_identity, same


def validate(queue, launch, result, allocation, recipe, recipe_sha, job, timezone):
    if type(job) is not int or job <= 0 or allocation["JobId"] != str(job):
        raise ValueError("Invalid or inconsistent job identity")
    source = [row for row in recipe["records"] if row["kind"] == "file" and
              row["path"] == str(RECIPE_ROOT / "benchmark_tools/submit_root_context_native_session.py")]
    if len(source) != 1 or any(not same(result["source"][key], source[0][key]) for key in ("sha256", "bytes")):
        raise ValueError("Waiting-session source differs from recipe")
    if not same(queue["command"], QUEUE_COMMAND) or type(queue["returncode"]) is not int or queue["returncode"] != 0 or queue["stdout"].strip():
        raise ValueError("Queue not verified empty")
    expected = dict(command=command(recipe_sha), recipe_sha256=recipe_sha, local_timeout_s=LOCAL_TIMEOUT,
        remote_timeout_s=REMOTE_TIMEOUT, remote_kill_grace_s=KILL_GRACE, scientific_timings_admitted=False)
    if any(not same(launch[k], v) for k, v in expected.items()):
        raise ValueError("Waiting command or bounds differ")
    exit_code, signal = map(int, allocation["ExitCode"].split(":"))
    expected = dict(status="wait_returned", returncode=1 if signal else exit_code, job_id=job,
        recipe_sha256=recipe_sha, scientific_timings_admitted=False, scheduler_terminal_verified=False)
    if any(not same(result[k], v) for k, v in expected.items()) or job_id(result["stdout"]) != job:
        raise ValueError("Waiting result differs from scheduler")
    times = [queue["started_unix_ns"], queue["finished_unix_ns"], launch["started_unix_ns"], result["finished_unix_ns"]]
    if any(type(t) is not int or t <= 0 for t in times) or times != sorted(times):
        raise ValueError("Receipt chronology differs")
    zone = ZoneInfo(timezone)
    parsed = [datetime.fromisoformat(allocation[key]) for key in ("StartTime", "EndTime")]
    if any(value.tzinfo is not None or value.microsecond for value in parsed):
        raise ValueError("Require second-resolution local scheduler timestamps")
    start, end = [int(value.replace(tzinfo=zone).timestamp()) * 10**9 for value in parsed]
    # Printed scheduler times have only second resolution; do not assert subsecond enclosure.
    if not (times[2] < start + 10**9 <= end + 10**9 and end < times[3]
            and times[3] - times[2] <= LOCAL_TIMEOUT * 10**9):
        raise ValueError("Waiting session does not enclose the job within its bound")
    return dict(status="native_bounded_session_verified", job_id=job, scheduler_timezone=timezone,
        printed_job_started_unix_ns=start, printed_job_ended_unix_ns=end,
        wait_started_unix_ns=times[2], wait_finished_unix_ns=times[3], scientific_timings_admitted=False,
        limitations=["Scheduler times have second resolution; controller/client wall clocks assumed consistent.",
            "Receipt binding does not validate native work, user-manager lifecycle, overhead or isolation."])


def audit(directory, scheduler_path, recipe_path, recipe_sha, job, timezone):
    paths = [directory / name for name in ("queue.json", "launch.json", "result.json")]
    evidence = [record(path) for path in [scheduler_path, recipe_path, *paths]]
    raw = scheduler_path.read_text()
    scheduler_identity(raw, job, require_completed=False)
    allocation = dict(re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", raw))
    recipe = read_pinned(recipe_path, recipe_sha)
    result = validate(*(json.loads(path.read_text()) for path in paths), allocation, recipe, recipe_sha, job, timezone)
    for item in evidence:
        check(item)
    return dict(result, evidence=evidence, source=record(__file__))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("directory", "scheduler", "recipe", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--recipe-sha", required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--scheduler-timezone", required=True)
    args = parser.parse_args()
    result = audit(args.directory.resolve(), args.scheduler.resolve(), args.recipe.resolve(),
                   args.recipe_sha, args.job, args.scheduler_timezone)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")

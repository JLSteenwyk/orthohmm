"""Verify a bounded waiting-session receipt against terminal Slurm evidence."""

import argparse
from datetime import datetime
import json
from pathlib import Path
from zoneinfo import ZoneInfo

from benchmark_tools.audit_root_context_controls import scheduler
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.submit_root_context_session import command, job_id, RECIPE


def validate(queue, launch, result, allocation, recipe, recipe_sha, job, timezone):
    source = [row for row in recipe["records"] if row["kind"] == "file"
              and row["path"] == str(RECIPE / "benchmark_tools/submit_root_context_session.py")]
    if (len(source) != 1 or result["source"]["sha256"] != source[0]["sha256"]
            or result["source"]["bytes"] != source[0]["bytes"]):
        raise ValueError("Waiting-session source differs from pinned recipe")
    if (queue["command"] != ["squeue", "-h", "-p", "spark", "-o", "%i %T %j"]
            or type(queue["returncode"]) is not int or queue["returncode"] != 0 or queue["stdout"].strip()):
        raise ValueError("Queue not verified empty")
    if (launch["command"] != command(recipe_sha) or launch["recipe_sha256"] != recipe_sha
            or launch["local_timeout_s"] != 1050 or launch["remote_timeout_s"] != 1020
            or launch["remote_kill_grace_s"] != 10 or launch["scientific_timings_admitted"] is not False):
        raise ValueError("Waiting command or bounds differ")
    exit_code, signal = map(int, allocation["ExitCode"].split(":"))
    expected_exit = 1 if signal else exit_code
    if (result["status"] != "wait_returned" or type(result["returncode"]) is not int
            or result["returncode"] != expected_exit or result["job_id"] != job or job_id(result["stdout"]) != job
            or result["recipe_sha256"] != recipe_sha or result["scientific_timings_admitted"] is not False
            or result["scheduler_terminal_verified"] is not False):
        raise ValueError("Waiting receipt disagrees with terminal scheduler")
    times = [queue["started_unix_ns"], queue["finished_unix_ns"], launch["started_unix_ns"], result["finished_unix_ns"]]
    if any(type(t) is not int or t <= 0 for t in times) or times != sorted(times):
        raise ValueError("Receipt chronology differs")
    zone = ZoneInfo(timezone)
    start, end = [int(datetime.fromisoformat(allocation[key]).replace(tzinfo=zone).timestamp()) * 10**9
                  for key in ("StartTime", "EndTime")]
    # Slurm's printed timestamps have one-second resolution.
    if not (times[2] < start + 10**9 <= end + 10**9 and end < times[3] and times[3]-times[2] <= 1050*10**9):
        raise ValueError("Waiting session does not enclose the scheduler job")
    return dict(status="bounded_session_submission_verified", job_id=job, scheduler_timezone=timezone,
        printed_job_started_unix_ns=start, printed_job_ended_unix_ns=end,
        wait_started_unix_ns=times[2], wait_finished_unix_ns=times[3], scientific_timings_admitted=False,
        limitations=["One-second scheduler timestamp resolution; local/controller wall clocks assumed consistent.",
            "Submission receipt and job lifetime only; workload, runtime, manager lifecycle and context replay are separate."])


def audit(directory, scheduler_path, recipe_path, recipe_sha, job, timezone):
    paths = [directory / name for name in ("queue.json", "launch.json", "result.json")]
    evidence = [record(path) for path in [*paths, scheduler_path, recipe_path]]
    values = [json.loads(path.read_text()) for path in paths]
    allocation = scheduler(scheduler_path.read_text(), job, "v2")
    recipe = read_pinned(recipe_path, recipe_sha)
    result = validate(*values, allocation, recipe, recipe_sha, job, timezone)
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

"""Validate recorded local allocation policy, not freshness or timing eligibility."""

import re
import math
import subprocess
import time

from benchmark_tools.probe_dgx_step_separation import save

from benchmark_tools.capture_job_scheduler import terminal_record


def validate(raw, job, phase, *, command, cwd, time_limit="1-00:00:00"):
    if type(job) is not int or job <= 0 or phase not in {"running", "terminal"}:
        raise ValueError("Require positive job ID and observation phase")
    if not all(isinstance(v, str) and v for v in (command, cwd, time_limit)):
        raise ValueError("Require explicit command, cwd and time limit")
    lines = [line for line in raw.splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError("Require one controller record")
    pairs = re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", lines[0])
    fields = dict(pairs)
    if len(fields) != len(pairs):
        raise ValueError("Duplicate scheduler field")
    # Exclusive-node allocation reserves all 192 slots, unlike the 64-slot task.
    expected = dict(JobId=str(job), Partition="gpu", NodeList="bizon",
        NumNodes="1", NumCPUs="192", NumTasks="1", OverSubscribe="NO",
        MinMemoryNode="128G", Requeue="0", Restarts="0", Command=command,
        WorkDir=cwd, TimeLimit=time_limit)
    expected["CPUs/Task"] = "64"
    if any(fields.get(k) != v for k, v in expected.items()) or any(
            k in fields for k in ("ArrayJobId", "ArrayTaskId", "HetJobId", "HetJobOffset")):
        raise ValueError("Threadripper allocation or execution identity differs")
    if not re.fullmatch(r"[0-9]+:[0-9]+", fields.get("ExitCode", "")):
        raise ValueError("Invalid scheduler exit status")
    if phase == "running":
        if fields.get("JobState") != "RUNNING":
            raise ValueError("Require running allocation")
    elif terminal_record(raw, job) is None:
        raise ValueError("Require terminal controller evidence")
    return dict(status="threadripper_controller_record_verified", job_id=job,
        scheduler_state=fields["JobState"], scheduler_exit_code=fields["ExitCode"],
        scheduler_terminal_verified=phase == "terminal", fields=fields,
        scientific_execution_authorized=False, scientific_timings_admitted=False,
        limitations=["Recorded controller fields only; caller must verify freshness and submission provenance.",
                     "Native affinity, effective cgroup limits and whole-host quietness require independent checks.",
                     "Caller must pin command/cwd/time limit in the production recipe; diagnostic records do not authorize production."])


def remaining_budget(raw, job, *, command, cwd, query_elapsed_s):
    allocation = validate(raw, job, "running", command=command, cwd=cwd,
                          time_limit="1-02:00:00")
    if (type(query_elapsed_s) not in (int, float) or not math.isfinite(query_elapsed_s)
            or not 0 <= query_elapsed_s <= 5):
        raise ValueError("Scheduler observation exceeded freshness allowance")
    match = re.fullmatch(r"(?:(\d+)-)?(\d{2}):([0-5]\d):([0-5]\d)",
                         allocation["fields"].get("RunTime", ""))
    if match is None or int(match[2]) > 23:
        raise ValueError("Invalid scheduler elapsed time")
    elapsed = int(match[1] or 0)*86400 + int(match[2])*3600 + int(match[3])*60 + int(match[4])
    # Allow for integer-second reporting, query latency and release scheduling.
    margin = 1 + math.ceil(query_elapsed_s) + 30
    available = 93600 - elapsed - margin
    if available < 85800 + 4200:
        raise ValueError("Insufficient allocation time for native timeout and reporting reserve")
    return dict(status="threadripper_release_budget_checked", allocation=allocation,
        scheduler_elapsed_s=elapsed, conservative_available_s=available,
        query_elapsed_s=query_elapsed_s, rounding_query_release_margin_s=margin,
        native_timeout_s=85800, post_native_reserve_s=4200,
        scientific_execution_authorized=False, scientific_timings_admitted=False)


class ReleaseBudgetGuard:
    """Fresh local query at the release gate; recipe and quiet-host checks are separate."""

    def __init__(self, job, *, command, cwd, runner=subprocess.run, clock=time.monotonic):
        self.job, self.command, self.cwd = job, command, cwd
        self.runner, self.clock = runner, clock

    def __call__(self, directory):
        argv = ["scontrol", "show", "job", str(self.job), "--oneliner"]
        observation = dict(command=argv, status="release_budget_check_failed")
        started = self.clock()
        try:
            result = self.runner(argv, capture_output=True, text=True, timeout=5)
            elapsed = self.clock() - started
            observation.update(returncode=result.returncode, stdout=result.stdout,
                               stderr=result.stderr, query_elapsed_s=elapsed)
            if result.returncode != 0:
                raise ValueError("Scheduler query failed")
            budget = remaining_budget(result.stdout, self.job, command=self.command,
                                      cwd=self.cwd, query_elapsed_s=elapsed)
            observation.update(status="release_budget_check_passed", budget=budget)
        except Exception as error:
            observation.update(error_type=type(error).__name__, error=str(error))
            save(directory / "release_budget.json", observation)
            raise
        save(directory / "release_budget.json", observation)
        return budget

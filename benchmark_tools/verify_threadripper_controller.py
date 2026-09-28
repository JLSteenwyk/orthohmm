"""Validate recorded local allocation policy, not freshness or timing eligibility."""

import re

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

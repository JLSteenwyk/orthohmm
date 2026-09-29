"""Read terminal Slurm accounting without treating missing usage as measured zero."""

import argparse
import csv
from datetime import datetime, timezone
import io
import json
from pathlib import Path
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import record

FIELDS = "JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,TotalCPU,MaxRSS,MaxVMSize,Start,End"


def assess(text, jobs, gather_type):
    reader = csv.DictReader(io.StringIO(text), delimiter="|")
    if reader.fieldnames != FIELDS.split(","):
        raise ValueError("Unexpected accounting columns")
    rows = list(reader)
    ids = [r["JobIDRaw"] for r in rows]
    if len(ids) != len(set(ids)) or not rows:
        raise ValueError("Missing or duplicate accounting rows")
    if any(j.split(".")[0] not in {str(n) for n in jobs} for j in ids):
        raise ValueError("Unexpected job")
    if any(str(j) not in ids for j in jobs):
        raise ValueError("Missing allocation row")
    disabled = gather_type in {None, "(null)", "jobacct_gather/none"}
    for row in rows:
        if None in row or any(v is None for v in row.values()):
            raise ValueError("Malformed accounting row")
        row["memory_usage_available"] = bool(row["MaxRSS"].strip()) and not disabled
        row["cpu_usage_available"] = bool(row["TotalCPU"].strip()) and not disabled
    return dict(rows=rows, jobacct_gather_type=gather_type,
        accounting_collection_configured=not disabled,
        missing_memory_rows=sum(not r["memory_usage_available"] for r in rows),
        unavailable_cpu_rows=sum(not r["cpu_usage_available"] for r in rows),
        complete_job_resource_accounting_verified=False, scientific_timings_admitted=False,
        limitations=["Current controller configuration is not proof of historical job-time configuration.",
            "Nonempty fields and an enabled collector do not establish accounting correctness or full-job scope.",
            "Step RSS maxima must not be summed or substituted for a simultaneous whole-job cgroup peak.",
            "Terminal state and elapsed time do not turn unavailable usage into measured zero."])


def capture(jobs, output):
    if output.exists():
        raise FileExistsError(output)
    if not jobs or len(jobs) != len(set(jobs)) or any(j <= 0 for j in jobs):
        raise ValueError("Require distinct positive job IDs")
    commands = [["sacct", "-j", ",".join(map(str, jobs)), "--parsable2", "--units=K", "--format=" + FIELDS],
                ["scontrol", "show", "config"]]
    evidence = []
    for command in commands:
        started = datetime.now(timezone.utc).isoformat()
        value = subprocess.run(command, capture_output=True, text=True, timeout=30, check=True)
        evidence.append(dict(command=command, started_utc=started, returncode=value.returncode,
                             stdout=value.stdout, stderr=value.stderr))
    config = dict(line.split("=", 1) for line in evidence[1]["stdout"].splitlines() if "=" in line)
    config = {k.strip(): v.strip() for k, v in config.items()}
    result = assess(evidence[0]["stdout"], jobs, config.get("JobAcctGatherType"))
    result.update(commands=evidence, source=record(Path(__file__).resolve()),
                  status="terminal_accounting_observation", publication_ready=False)
    with output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", type=int, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = capture(args.job, args.output)
    print(json.dumps({k: result[k] for k in ("missing_memory_rows", "unavailable_cpu_rows", "jobacct_gather_type")}))

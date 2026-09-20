"""Check replacement-scaling controller records without authorizing execution."""

import argparse
import json
from pathlib import Path
import re

from benchmark_tools.capture_job_scheduler import terminal_record
from benchmark_tools.measure_root_context_scaling import load_task, PLAN_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.verify_scaling_task_records import RECIPE_ROOT

SUBMISSION_SCRIPT = RECIPE_ROOT / "benchmark_tools/run_dgx_root_context_scaling.sh"


def validate(raw, job, phase):
    if type(job) is not int or job <= 0:
        raise ValueError("Require positive integer job identity")
    if phase not in ("running", "terminal"):
        raise ValueError("Require running or terminal observation phase")
    lines = [line for line in raw.splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError("Require exactly one oneline controller record")
    pairs = re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", lines[0])
    fields = dict(pairs)
    if len(fields) != len(pairs):
        raise ValueError("Duplicate scheduler field")
    expected = dict(JobId=str(job), Partition="spark", Restarts="0", Requeue="0",
        NodeList="spark-7ff0", OverSubscribe="NO", MinMemoryNode="96G",
        NumNodes="1", NumCPUs="20", NumTasks="1", WorkDir=str(RECIPE_ROOT),
        Command=str(SUBMISSION_SCRIPT))
    expected["CPUs/Task"] = "20"
    if any(fields.get(k) != v for k, v in expected.items()) or any(
            k in fields for k in ("ArrayJobId", "ArrayTaskId", "HetJobId", "HetJobOffset")):
        raise ValueError("Allocation identity or matched resource policy differs")
    # Slurm can print a day as either an hours count or a day-prefixed duration.
    if fields.get("TimeLimit") not in ("1-00:00:00", "24:00:00"):
        raise ValueError("Require the frozen 24-hour allocation limit")
    if not re.fullmatch(r"[0-9]+:[0-9]+", fields.get("ExitCode", "")):
        raise ValueError("Missing or invalid scheduler exit status")
    if phase == "running":
        if fields.get("JobState") != "RUNNING":
            raise ValueError("Require running allocation; pending/completing is insufficient")
    elif terminal_record(raw, job) is None:
        raise ValueError("Require terminal controller evidence; timeout of observation is insufficient")
    return dict(status="scaling_allocation_record_verified", job_id=job, phase=phase,
        fields=fields, scheduler_terminal_verified=phase == "terminal",
        scientific_timings_admitted=False, execution_authorized=False,
        next_submission_authorized=False, automatic_retry=False,
        limitations=["Checks recorded allocation identity and resources, not record freshness or authenticity.",
            "Does not bind a task index, submitted arguments, source recipe, SSH session or native result.",
            "Exclusive Slurm allocation does not exclude services or work outside Slurm.",
            "Environmental policy, whole-run observations and output/failure audits remain required."])


def audit(plan_path, index, scheduler_path, job, phase):
    evidence = [record(plan_path), record(scheduler_path)]
    plan, task, _ = load_task(plan_path, index)
    expected = dict(cpus=20, exclusive=True, max_concurrent_runs=1, memory_gib=96,
        node="spark-7ff0", partition="spark", requeue=False, time_limit_s=86400)
    if plan["allocation"] != expected:
        raise ValueError("Frozen allocation specification differs")
    result = validate(Path(scheduler_path).read_text(), job, phase)
    for item in evidence:
        check(item)
    return dict(result, plan_sha256=PLAN_SHA, requested_index=task["index"],
        task_identity_bound=False, evidence=evidence, source=record(__file__))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "scheduler", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--phase", choices=("running", "terminal"), required=True)
    args = parser.parse_args(argv)
    result = audit(args.plan, args.index, args.scheduler, args.job, args.phase)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

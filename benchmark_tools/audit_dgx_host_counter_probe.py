"""Validate scheduled counter availability without admitting controlled timing."""

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.probe_host_counters import summarize
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.slurm_resource_snapshot import scoped_path, counters

PROBE_SHA = "3d3c193305085e17aa10e05be842d814f66d020d7dde285c8d9288ead2ce3f3c"


def validate(data, scheduler, job):
    expected = {"JobId": str(job), "JobState": "COMPLETED", "ExitCode": "0:0",
        "NodeList": "spark-7ff0", "Partition": "spark", "CPUs/Task": "1",
        "NumCPUs": "20", "MinMemoryNode": "1G", "OverSubscribe": "NO", "Restarts": "0"}
    fields = {}
    for key, value in expected.items():
        matches = re.findall(r"(?:^|\s)" + re.escape(key) + r"=([^\s]+)", scheduler)
        if matches != [value]:
            raise ValueError("Unexpected scheduler field: " + key)
        fields[key] = value
    if (data["status"] != "prospective_host_counter_probe" or data["hostname"] != "spark-7ff0"
            or data["source"]["sha256"] != PROBE_SHA or data["publication_ready"] is not False
            or len(data["snapshots"]) != 2):
        raise ValueError("Wrong probe identity")
    for sample in data["snapshots"]:
        scope = scoped_path(sample["raw"]["cgroup_membership"], job)
        if "step_batch" not in scope.parts or sample["errors"]:
            raise ValueError("Wrong scope or incomplete counter availability")
        required = {"cgroup_cpu.stat", "cgroup_memory.current", "cgroup_memory.peak", "cgroup_memory.events",
                    "host_cpu_pressure", "host_memory_pressure", "host_io_pressure"}
        if set(sample["optional"]) != required:
            raise ValueError("Missing counter fields")
        if not {"usage_usec", "user_usec", "system_usec"} <= counters(sample["optional"]["cgroup_cpu.stat"]).keys():
            raise ValueError("Missing cumulative cgroup CPU accounting")
    replay = summarize(*data["snapshots"], data["clock_ticks_per_second"])
    if replay != data["summary"]:
        raise ValueError("Raw replay differs")
    return fields, str(scope)


def run(probe, job, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    probe_record = record(probe)
    data = json.loads(probe.read_text())
    command = ["scontrol", "show", "job", str(job)]
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    fields, scope = validate(data, result.stdout, job)
    check(probe_record)
    report = dict(status="scheduled_host_counter_availability_verified", publication_ready=False,
        controlled_workload_verified=False, job_id=job, probe=probe_record, source=record(__file__),
        scheduler_command=command, scheduler_raw=result.stdout, scheduler_stderr=result.stderr,
        scheduler_fields=fields, scope=scope, summary=data["summary"],
        limitations=["Read-only batch-step availability check; not native-process measurement or overhead calibration.",
            "Exclusive allocation reserves 20 CPUs; requested CPUs per task is 1. No load was generated.",
            "Does not validate foreign-CPU estimation, short-lived-load detection or past timing admission.",
            "Scheduler times are retained unqualified; no timezone inference is made."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.probe.resolve(), args.job, args.output.absolute())

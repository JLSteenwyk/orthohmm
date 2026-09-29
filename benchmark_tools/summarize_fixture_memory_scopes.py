"""Expose retained diagnostic memory scopes without granting timing admission."""

import argparse
import json
from pathlib import Path

from benchmark_tools.build_private_timing_environment import record, write
from benchmark_tools.prepare_ob_candidate_neighborhood import check
from benchmark_tools.slurm_resource_snapshot import counters, scoped_path


def native_cpu(done, job_id, completion):
    """Report raw counter differences, including the wrapper's bracket work."""
    before, after = done["snapshots"]
    times = [before["started_monotonic_ns"], before["finished_monotonic_ns"],
             done["started_ns"], done["finished_ns"],
             after["started_monotonic_ns"], after["finished_monotonic_ns"]]
    if (any(type(t) is not int or t <= 0 for t in times)
            or not times[0] <= times[1] < times[2] < times[3] < times[4] <= times[5]):
        raise ValueError("CPU reads must bracket command execution")
    if completion["status"] != "anchor_only_at_boundaries" or completion["errors"]:
        raise ValueError("Unverified native completion")
    if completion["command_finished_ns"] != done["finished_ns"]:
        raise ValueError("Completion and command disagree")
    for key in ("boot_id", "online_cpus", "cgroup_membership"):
        if not before["raw"][key] or before["raw"][key] != after["raw"][key]:
            raise ValueError("CPU scope or host changed")
    scope = scoped_path(before["raw"]["cgroup_membership"], job_id)
    user_path = Path(completion["before"]["scope"])
    user_scope = Path("/") / user_path.relative_to("/sys/fs/cgroup")
    if (not scope.is_relative_to(user_scope)
            or completion["after"]["scope"] != str(user_path)):
        raise ValueError("CPU scope is outside native completion subtree")
    values = []
    fields = ("usage_usec", "user_usec", "system_usec")
    for sample in (before, after):
        if sample["errors"]:
            raise ValueError("Incomplete native CPU observation")
        parsed = counters(sample["optional"]["cgroup_cpu.stat"])
        if not set(fields) <= parsed.keys():
            raise ValueError("Missing native CPU fields")
        values.append({key: parsed[key] for key in fields})
    delta = {key: values[1][key] - values[0][key] for key in fields}
    if any(value < 0 for value in delta.values()):
        raise ValueError("Native CPU counter decreased")
    return dict(scope=str(scope), before=values[0], after=values[1], delta_usec=delta,
                cpu_seconds=delta["usage_usec"] / 1e6,
                command_wall_seconds=(times[3] - times[2]) / 1e9,
                read_window_ns=times, final_whole_job_cpu_seconds=None,
                scientific_timings_admitted=False,
                limitations=["CPU differences include wrapper work within the read bracket.",
                             "Not terminal whole-job CPU, pure algorithm CPU or controlled timing.",
                             "Completion observations do not prove continuous containment."])


def peak(observation):
    value = observation["raw"]["memory.peak"].strip()
    if observation["errors"] or not value.isascii() or not value.isdecimal():
        raise ValueError("Incomplete or invalid memory peak")
    if (type(observation["started_ns"]) is not int or type(observation["finished_ns"]) is not int
            or not 0 < observation["started_ns"] <= observation["finished_ns"]):
        raise ValueError("Invalid read interval")
    return dict(bytes=int(value), scope=observation["scope"],
                read_started_ns=observation["started_ns"], read_finished_ns=observation["finished_ns"])


def summarize(replayed):
    if (replayed["status"] != "threadripper_scaling_measurement_replayed"
            or replayed["native_completion"]["status"] != "anchor_only_at_boundaries"):
        raise ValueError("Require completed v5 replay")
    job = replayed["job_memory"]
    values = dict(native_step=peak(replayed["memory"]),
                  job_before_native=peak(job["before"]), job_after_native=peak(job["after"]),
                  job_through_reporting=peak(replayed["report_finalization"]["job_memory"]))
    parent = values["job_before_native"]["scope"]
    if (any(values[k]["scope"] != parent for k in ("job_after_native", "job_through_reporting"))
            or Path(values["native_step"]["scope"]).parent != Path(parent)):
        raise ValueError("Memory scopes are not the expected parent and step")
    if not (values["job_before_native"]["bytes"] <= values["job_after_native"]["bytes"]
            <= values["job_through_reporting"]["bytes"]
            and values["native_step"]["bytes"] <= values["job_after_native"]["bytes"]):
        raise ValueError("Memory peaks contradict scope or monotonicity")
    return dict(measurements=values, final_whole_job_peak_bytes=None,
                scientific_timings_admitted=False,
                limitations=["Native step includes its launcher; not pure algorithm RSS.",
                             "Job peaks include preparation and observation only through their read boundaries.",
                             "Overlapping peaks must not be added or subtracted.",
                             "Final whole-job peak is not established by these receipts."])


def run(repo, output):
    rows = []
    for job_id in (22373, 22374, 22375):
        path = repo / f"benchmark_tools/results/threadripper_private_collector_{job_id}.json"
        receipt = json.loads(path.read_text())
        ref = receipt.get("replay", receipt.get("audit"))
        check(ref)
        report = json.loads(Path(ref["path"]).read_text())
        replayed = report.get("replay", report)
        for evidence in replayed["evidence"]:
            check(evidence)
        row = summarize(replayed)
        done_refs = [item for item in replayed["evidence"] if Path(item["path"]).name == "done.json"]
        if len(done_refs) != 1:
            raise ValueError("Require one retained native command receipt")
        done = json.loads(Path(done_refs[0]["path"]).read_text())
        row["native_cpu"] = native_cpu(done, job_id, replayed["native_completion"])
        row["native_command_source"] = done_refs[0]
        for evidence in replayed["evidence"]:
            check(evidence)
        row.update(job_id=job_id, source=record(path), replay_source=ref)
        rows.append(row)
    result = dict(status="diagnostic_memory_scopes_summarized", rows=rows, source=record(Path(__file__)),
                  scientific_timings_admitted=False, production_runs_launched=0)
    write(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.output.resolve())["status"])

"""Reproduce retained DGX resource accounting without certifying scientific timing."""

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.gnu_time_companion import parse as parse_time
from benchmark_tools.measure_slurm_command import check_allocation
from benchmark_tools.monitor_slurm_resources import summarize
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.slurm_resource_snapshot import counters, cpus, scoped_path


def validate_rows(measured, rows):
    if len(rows) < 2 or measured["status"] != "command_exited_zero" or measured["timed_out"] is not False:
        raise ValueError("Require complete measurement and bracketing samples")
    domain = measured["clock_domain"]
    if (domain["hostname"] != "spark-7ff0" or domain["clock"] != "time.monotonic"
            or domain["unit"] != "seconds" or not domain["boot_id"]):
        raise ValueError("Unknown resource clock domain")
    start, launch, end, finish = (measured[k] for k in ("wrapper_started_monotonic_s",
        "command_launch_started_monotonic_s", "command_wait_finished_monotonic_s", "wrapper_finished_monotonic_s"))
    if not all(math.isfinite(v) for v in (start, launch, end, finish)) or not 0 < start <= launch < end <= finish:
        raise ValueError("Invalid command/wrapper clock boundaries")
    if measured["command_wall_s"] != end - launch or measured["wrapper_wall_s"] != finish - start:
        raise ValueError("Recorded wall durations differ from clock boundaries")
    previous_finish = start
    for index, row in enumerate(rows):
        begin, done = row["started_monotonic_s"], row["finished_monotonic_s"]
        if (row["index"] != index or not all(math.isfinite(v) for v in (begin, done, row["elapsed_s"]))
                or not previous_finish <= begin <= done <= finish or row["elapsed_s"] != done - start):
            raise ValueError("Invalid observation chronology")
        previous_finish = done
        sample = row["snapshot"]
        if (sample["job_id"] != measured["job_id"] or sample["pid"] != measured["wrapper_pid"]
                or sample["scope"] != measured["baseline_scope"] or sample["status"] != "resource_snapshot"
                or str(scoped_path(sample["proc_cgroup"], measured["job_id"])) != sample["scope"]):
            raise ValueError("Resource snapshot identity differs")
        metrics, raw = sample["metrics"], sample["metrics"]["raw"]
        expected = {"cpu": counters(raw["cpu.stat"]), "memory_events": counters(raw["memory.events"]),
            "memory_stat": counters(raw["memory.stat"]), "effective_cpus": cpus(raw["cpuset.cpus.effective"]),
            "memory_current_bytes": int(raw["memory.current"]),
            "memory_peak_since_creation_or_reset_bytes": int(raw["memory.peak"]),
            "memory_max_bytes": None if raw["memory.max"].strip() == "max" else int(raw["memory.max"])}
        if any(metrics[k] != v for k, v in expected.items()):
            raise ValueError("Parsed resource metrics differ from raw counters")
        processes = sample["processes"]
        if len({p["pid"] for p in processes}) != len(processes):
            raise ValueError("Duplicate process RSS observation")
        if (any(type(p["rss_bytes"]) is not int or p["rss_bytes"] < 0
                or not Path(p["cgroup"]).is_relative_to(Path(sample["scope"])) for p in processes)
                or sum(p["rss_bytes"] for p in processes) != sample["sampled_sum_process_rss_bytes"]):
            raise ValueError("Recorded process RSS sum or ownership differs")
        if len(metrics["effective_cpus"]) != measured["requested_cpus"]:
            raise ValueError("Resource cpuset changed")
    if rows[0]["finished_monotonic_s"] > launch or rows[-1]["started_monotonic_s"] < end:
        raise ValueError("Resource samples do not bracket command")
    check_allocation(rows[0]["snapshot"], measured["wrapper_pid"],
                     measured["requested_cpus"], measured["requested_memory_bytes"])
    if any(p["pid"] != measured["wrapper_pid"] for p in rows[-1]["snapshot"]["processes"]):
        raise ValueError("Native processes remain at final observation")
    summary = summarize(rows)
    if summary != measured["summary"]:
        raise ValueError("Resource summary does not reproduce")
    return summary


def replay(directory):
    paths = [directory / "measurement/results.json", directory / "measurement/samples.jsonl", directory / "native.time.tsv"]
    evidence = [record(p) for p in paths]
    measured = json.loads(paths[0].read_text())
    if measured["samples.jsonl"]["sha256"] != evidence[1]["sha256"]:
        raise ValueError("Raw resource series hash differs")
    helpers = [record(Path(__file__).with_name(name)) for name in
               ("measure_slurm_command.py", "slurm_resource_snapshot.py", "monitor_slurm_resources.py", "gnu_time_companion.py")]
    if measured["source"]["sha256"] != helpers[0]["sha256"] or measured["snapshot_source"]["sha256"] != helpers[1]["sha256"]:
        raise ValueError("Deployed resource source differs from replay source")
    with paths[1].open() as stream:
        rows = [json.loads(line) for line in stream]
    summary = validate_rows(measured, rows)
    timing = parse_time(paths[2].read_text())
    if timing["exit_status"] != measured["exit_code"]:
        raise ValueError("GNU-time exit status differs")
    cpu_s = summary["cpu_delta_usec"]["usage_usec"] / 1e6
    native_cpu_s = timing["user_seconds"] + timing["system_seconds"]
    for item in [*evidence, *helpers]:
        check(item)
    return {"status": "resource_accounting_reproduced_not_timing_admission", "evidence": evidence,
        "helpers": helpers, "source": record(__file__), "job_id": measured["job_id"],
        "clock_domain": measured["clock_domain"], "summary": summary, "gnu_time": timing,
        "reconciliation": {"cgroup_observed_cpu_seconds": cpu_s, "gnu_time_native_cpu_seconds": native_cpu_s,
            "cgroup_minus_gnu_cpu_seconds": cpu_s - native_cpu_s,
            "collector_command_wall_seconds": measured["command_wall_s"],
            "collector_minus_gnu_wall_seconds": measured["command_wall_s"] - timing["elapsed_seconds"],
            "maximum_sample_gap_seconds": max(b["started_monotonic_s"] - a["finished_monotonic_s"] for a, b in zip(rows, rows[1:])),
            "baseline_cgroup_peak_bytes": rows[0]["snapshot"]["metrics"]["memory_peak_since_creation_or_reset_bytes"]},
        "scientific_timing_admitted": False, "controlled_workload_verified": False,
        "limitations": ["Cgroup deltas include collector overhead and pre/post command observation margins.",
            "GNU-time CPU covers native child and waited-for descendants, with its own rounding and scope.",
            "Differences are reported without a fitted acceptance threshold or correction.",
            "Sampled aggregate RSS, cgroup memory peak and maximum process RSS are distinct statistics.",
            "Resource reproduction does not validate host isolation or native biological outputs."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    results = []
    for index in range(27):
        try:
            results.append({"index": index, **replay(args.root / f"run_{index:02d}")})
        except Exception as error:
            results.append({"index": index, "status": "replay_failed", "error_type": type(error).__name__, "error": str(error)})
    report = {"source": record(__file__), "runs": results, "failures": sum(r["status"] == "replay_failed" for r in results),
              "scientific_timings_admitted": 0}
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"runs": len(results), "failures": report["failures"]}))

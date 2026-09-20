"""Describe retained memory/I/O observations without certifying isolation."""

import argparse
import gzip
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.slurm_resource_snapshot import counters

AUDIT_SHA = "6da8441923e51a404fb2e4fec4c7302593fc1a1bc7e338b8996b9d2d2f24641a"


def describe(audit):
    runs = audit["runs"]
    if (type(audit["validated_tasks"]) is not int or audit["validated_tasks"] != 18
            or audit["issues"] != [] or len(runs) != 18
            or audit["scientific_timings_admitted"] is not False):
        raise ValueError("Require the complete audited panel without timing admission")
    rows = []
    for index, run in enumerate(runs):
        if (type(run["index"]) is not int or run["index"] != index or run["status"] != "validated"
                or run["scientific_timings_admitted"] is not False):
            raise ValueError("Incomplete, reordered or invalid outcomes")
        pressure = run["pressure_observation_window"]
        if (pressure["status"] != "native_step_pressure_diagnostic"
                or pressure["scientific_timings_admitted"] is not False
                or pressure["controlled_workload_verified"] is not False):
            raise ValueError("Unexpected pressure provenance or isolation claim")
        totals = pressure["native_stall_usec"]
        if set(totals) != {"cpu", "memory", "io"}:
            raise ValueError("Incomplete pressure resources")
        for values in totals.values():
            if (set(values) != {"some", "full"}
                    or any(type(v) is not int or v < 0 for v in values.values())):
                raise ValueError("Invalid stall counters")
        memory = run["memory"]
        if memory["errors"]:
            raise ValueError("Memory observation errors")
        events = counters(memory["raw"]["memory.events"])
        if not {"low", "high", "max", "oom", "oom_kill"} <= events.keys():
            raise ValueError("Incomplete memory events")
        values = [memory["raw"][key].strip() for key in ("memory.current", "memory.peak")]
        if any(not v.isascii() or not v.isdecimal() for v in values):
            raise ValueError("Invalid memory gauge")
        current, peak = map(int, values)
        if current > peak:
            raise ValueError("Current memory exceeds retained peak")
        rows.append(dict(index=index, method=run["method"], arm=run["arm"],
            native_stall_usec=totals, cgroup_memory_peak_bytes=peak,
            final_cgroup_memory_current_bytes=current, memory_events=events,
            original_flags=len(run["original_flagged_intervals"]),
            narrow_flags=len(run["narrow_flagged_intervals"])))
    return dict(status="retained_non_cpu_observations_described", runs=rows,
        observed_memory_stall_tasks=[r["index"] for r in rows if any(r["native_stall_usec"]["memory"].values())],
        observed_io_stall_tasks=[r["index"] for r in rows if any(r["native_stall_usec"]["io"].values())],
        nonzero_memory_event_tasks=[r["index"] for r in rows if any(r["memory_events"].values())],
        scientific_timings_admitted=False, environmental_validity_established=False,
        limitations=["Post-outcome description, not a prospective inclusion rule or a new native experiment.",
            "Native-step pressure covers the observation window and includes wrappers, not exact command-only stalls.",
            "Pressure totals cannot identify the cause of stalls or be subtracted from host pressure.",
            "Zero observed memory stalls/events do not establish memory-bandwidth or cache isolation.",
            "Cgroup peak includes charged cache/kernel memory; it is not maximum-process RSS.",
            "GPU activity, device-level I/O attribution and thermal interference are not established by these records.",
            "All original/narrow flags remain; no thresholds, timings or admission decisions change."])


def run(path):
    sources = [record(path), record(__file__), record(Path(__file__).with_name("slurm_resource_snapshot.py"))]
    if sources[0]["sha256"] != AUDIT_SHA:
        raise ValueError("Unexpected source audit hash")
    with gzip.open(path, "rt") as stream:
        result = describe(json.load(stream))
    for source in sources:
        check(source)
    return dict(result, sources=sources)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = run(args.audit.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")

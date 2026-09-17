"""Replay and summarize the complete frozen DGX overhead panel without exclusions."""

import argparse
import csv
import io
import json
from pathlib import Path
import re
import statistics
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_slurm_measurement import audit
from benchmark_tools.monitor_slurm_resources import source_record

MODES = ("sparse", "sampled", "sampled", "sparse", "sparse", "sampled")
HASHES = (
    "4e2d0066702b61fa51ef8612659ad6bef191cc937e74aca670d6d813fdc5fe06",
    "639428cc9057ecf7bf39327e5965030da5eb2a0334b11c02a1ef2c27b41f1f26",
    "dc68aca25d5b281175d682d5b4650d3c42f675bf0da78d4e493c80270a18d017",
    "5fff867b3281dd95b054f83f6cd1f8206c43e2423da83a4df6eb2d3d12427c0e",
    "78fd9d9b8f50989dd57687e625ace0065f12afbd7d918b35ff2b147a35742b89",
    "789bfaed9d5b5698bfa193980e49ee2a1afb301ba17fba589b85f65ee4fab764",
)
ROOT = "/home/jlsteenwyk/projects/orthohmm-publication/collector_load_recipe_v2"


def paired_summary(rows):
    if [r["mode"] for r in rows] != list(MODES) or [r["index"] for r in rows] != list(range(6)):
        raise ValueError("Changed complete counterbalanced inventory")
    pairs = []
    for pair, (sparse, sampled) in enumerate(((0, 1), (3, 2), (4, 5))):
        a, b = rows[sparse], rows[sampled]
        pairs.append({"pair": pair, "sparse_index": sparse, "sampled_index": sampled,
                      "wall_inflation": b["wall_s"] / a["wall_s"] - 1,
                      "cpu_difference_s": b["cpu_s"] - a["cpu_s"]})
    values = [p["wall_inflation"] for p in pairs]
    median = statistics.median(values)
    return {"pairs": pairs, "median_wall_inflation": median, "range_wall_inflation": [min(values), max(values)],
            "numerical_wall_budget_met": median <= 0.05 and max(values) <= 0.10}


def report(root, protocol):
    if source_record(protocol)["sha256"] != "eb28b990cc52dadf20ccf7d9201eb70abdaa86aa5c7b8220ec8d87c0ee0a650a":
        raise ValueError("Changed overhead protocol")
    accounting = subprocess.check_output(["sacct", "-j", "21640", "--parsable2", "--noheader",
        "--format=JobID,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)
    jobs = {}
    for row in csv.reader(io.StringIO(accounting), delimiter="|"):
        if re.fullmatch(r"21640_[0-5]", row[0]):
            if row[0] in jobs:
                raise ValueError("Duplicate scheduler identity")
            jobs[row[0]] = row
    if set(jobs) != {f"21640_{i}" for i in range(6)}:
        raise ValueError("Incomplete scheduler inventory")
    rows, checksum, previous_finish, domain = [], None, None, None
    for index, mode in enumerate(MODES):
        job = jobs[f"21640_{index}"]
        if job[1:3] != ["COMPLETED", "0:0"] or job[4] != "20" or job[5] not in {"96G", "96Gn"} or job[6] != "spark-7ff0":
            raise ValueError("Scheduler completion/allocation differs")
        directory = root / f"collector_load_panel_v1_{index}_{mode}"
        replay = audit(directory, HASHES[index])
        r = json.loads((directory / "results.json").read_text())
        if r["command"] != [ROOT + "/load", "8000000000", "20"] or r["cwd"] != ROOT:
            raise ValueError("Wrong workload command")
        expected_interval = 1 if mode == "sampled" else 86400
        if r["interval_s"] != expected_interval or r["host_interval_s"] != 30 or r["timeout_s"] != 600:
            raise ValueError("Changed observation cadence")
        if r["requested_cpus"] != 20 or r["requested_memory_bytes"] != 96 * 1024 ** 3:
            raise ValueError("Wrong resource allocation")
        if domain is None:
            domain = r["clock_domain"]
        if r["clock_domain"] != domain or (previous_finish is not None and previous_finish > r["wrapper_started_monotonic_s"]):
            raise ValueError("Changed clock domain or overlapping runs")
        previous_finish = r["wrapper_finished_monotonic_s"]
        native = json.loads((directory / "command.log").read_text())
        if native["workers"] != 20 or native["iterations_per_worker"] != 8000000000 or len(native["checksums"]) != 20:
            raise ValueError("Incorrect workload completion")
        if any(not re.fullmatch(r"[0-9a-f]{16}", value) for value in native["checksums"]):
            raise ValueError("Invalid native checksum")
        if checksum is None:
            checksum = native
        if native != checksum:
            raise ValueError("Workload checksums differ")
        host_rows = [json.loads(line) for line in (directory / "host_samples.jsonl").read_text().splitlines()]
        durations = [row["snapshot"]["finished_monotonic_s"] - row["snapshot"]["started_monotonic_s"] for row in host_rows]
        host = r["host_workload"]
        gates = {"minimum_duration": r["command_wall_s"] >= 60,
                 "full_load": r["summary"]["mean_cpu_cores_over_observed_span"] >= 18,
                 "host_cadence": len(host_rows) >= (3 if mode == "sampled" else 2),
                 "low_observed_competition": host["maximum_observed_foreign_average_cores"] < 0.25,
                 "no_host_errors": host["observation_errors"] == 0}
        rows.append({"index": index, "mode": mode, "scheduler": job, "replay": replay,
                     "wall_s": r["command_wall_s"], "cpu_s": r["summary"]["cpu_delta_usec"]["usage_usec"] / 1e6,
                     "mean_cpu_cores": r["summary"]["mean_cpu_cores_over_observed_span"],
                     "max_sampled_sum_rss_bytes": r["summary"]["maximum_sampled_sum_rss_bytes"],
                     "max_reported_cgroup_peak_bytes": r["summary"]["maximum_reported_cgroup_peak_bytes"],
                     "host_scan_durations_s": durations, "max_foreign_cores": host["maximum_observed_foreign_average_cores"],
                     "max_host_gap_s": host["maximum_between_snapshot_gap_s"], "gates": gates})
    result = paired_summary(rows)
    return {**result, "status": "bounded_overhead_panel_evaluated", "runs": rows, "native_work": checksum,
            "observed_load_and_workload_gates_met": all(all(r["gates"].values()) for r in rows),
            "protocol": source_record(protocol), "source": source_record(__file__), "accounting": accounting,
            "scientific_execution_authorized": False, "all_protocol_identity_gates_verified": False,
            "limitations": ["Sparse mode still has pre/post observations; this is incremental periodic-observation cost.",
                            "Three counterbalanced pairs are descriptive, not a significance test or an overhead correction.",
                            "One compute process,20 threads,small arenas: not disk-heavy or large-process-inventory behavior.",
                            "The original runner pins native/collector files but does not record full interpreter/environment identity before and after each task.",
                            "Complete scientific runtime freeze remains required; this report does not waive that identity gate.",
                            "No overhead subtraction or exclusion of unfavorable runs is authorized."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "protocol", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = report(args.root.resolve(), args.protocol.resolve())
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")

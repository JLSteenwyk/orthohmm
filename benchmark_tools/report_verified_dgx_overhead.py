"""Audit the complete prospective verified-wrapper DGX overhead repeat."""

import argparse
import csv
import io
import json
import math
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_slurm_measurement import audit
from benchmark_tools.gnu_time_companion import command, parse
from benchmark_tools.monitor_slurm_resources import source_record
from benchmark_tools.report_dgx_collector_overhead import MODES, ROOT, paired_summary

PROJECT = str(Path(ROOT).parent)
PROTOCOL_SHA = "d04f2c114fdf4934cb3cff20d511f0678f1abb3f2873adb436a636e1ca0cda87"
RUNTIMES = (
    ("trees.json", "2f38fc57683e51a7b6768b4709db16293590983640c41cfe12a56ed6036750ae", 26673),
    ("system_trees.json", "4083d15c0ffc40013756c4588763aa8af24ca39165622665ee5074662e4b46ce", 10066),
    ("recipe_trees_v2.json", "ee0f6b5befb41acb4627e5d9a3c00b1c9490879c7bb2075997bd1c933960fbb8", 5),
)


def completed_jobs(accounting):
    jobs = {}
    for row in csv.reader(io.StringIO(accounting), delimiter="|"):
        if row and re.fullmatch(r"21647_[0-5]", row[0]):
            if row[0] in jobs:
                raise ValueError("Duplicate scheduler identity")
            jobs[row[0]] = row
    if set(jobs) != {f"21647_{i}" for i in range(6)}:
        raise ValueError("Incomplete scheduler inventory")
    for row in jobs.values():
        if (row[1:3] != ["COMPLETED", "0:0"] or row[4] != "20"
                or row[5] not in {"96G", "96Gn"} or row[6] != "spark-7ff0"):
            raise ValueError("Panel is incomplete, failed or incorrectly allocated")
    return jobs


def check_verification(verified, measured):
    expected = [{"path": PROJECT + "/runtime_inventory_v1/" + name, "sha256": sha,
                 "records": count, "scientific_execution_authorized": False,
                 "status": "runtime_tree_identity_matches"} for name, sha, count in RUNTIMES]
    if (verified["status"] != "command_exited_zero" or verified["before"] != expected
            or verified["after"] != expected or verified["measurement"] != measured
            or verified["source_sha256"] != "36243fcdfe6f49e73dbe1dd1b627bf4a330f19f3cd52dcf94a9618ddf0415cf8"):
        raise ValueError("Before/after identity or embedded collector record differs")
    for key in ("before_check_wall_s", "after_check_wall_s"):
        if not math.isfinite(verified[key]) or verified[key] <= 0:
            raise ValueError("Invalid verification duration")


def report(root, protocol):
    protocol_record = source_record(protocol)
    if protocol_record["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed prospective protocol")
    accounting = subprocess.check_output(["sacct", "-j", "21647", "--parsable2", "--noheader",
        "--format=JobID,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)
    jobs = completed_jobs(accounting)  # Never inspect partial panel outcomes.
    rows, native_work, clock_domain, previous_end = [], None, None, None
    for index, mode in enumerate(MODES):
        name = f"verified_collector_panel_v2_{index}_{mode}"
        directory = root / name
        paths = [directory / "verification.json", directory / "measurement/results.json",
                 directory / "native.time.tsv"]
        records = [source_record(path) for path in paths]
        verified, measured = (json.loads(path.read_text()) for path in paths[:2])
        check_verification(verified, measured)
        replay = audit(directory / "measurement", records[1]["sha256"])
        expected = command([ROOT + "/load", "8000000000", "20"], PROJECT + "/" + name + "/native.time.tsv")
        if (measured["command"] != expected or measured["cwd"] != ROOT
                or measured["status"] != "command_exited_zero" or measured["exit_code"] != 0
                or measured["timed_out"]):
            raise ValueError("Incorrect command or unsuccessful collector completion")
        if (measured["interval_s"] != (1 if mode == "sampled" else 86400)
                or measured["host_interval_s"] != 30 or measured["timeout_s"] != 600
                or measured["requested_cpus"] != 20 or measured["requested_memory_bytes"] != 96 * 1024 ** 3):
            raise ValueError("Changed cadence or allocation")
        if clock_domain is None:
            clock_domain = measured["clock_domain"]
        if (measured["clock_domain"] != clock_domain or (previous_end is not None
                and previous_end > measured["wrapper_started_monotonic_s"])):
            raise ValueError("Changed clock domain or overlapping measurement spans")
        previous_end = measured["wrapper_finished_monotonic_s"]
        native = json.loads((directory / "measurement/command.log").read_text())
        if (native["workers"] != 20 or native["iterations_per_worker"] != 8000000000
                or len(native["checksums"]) != 20
                or any(not re.fullmatch(r"[0-9a-f]{16}", value) for value in native["checksums"])):
            raise ValueError("Incomplete native workload")
        if native_work is None:
            native_work = native
        if native != native_work:
            raise ValueError("Native workload checksums differ")
        timing = parse(paths[2].read_text())
        if timing["exit_status"] != 0:
            raise ValueError("GNU-time command failed")
        host = measured["host_workload"]
        summary = measured["summary"]
        gates = {"duration": measured["command_wall_s"] >= 60,
                 "full_load": summary["mean_cpu_cores_over_observed_span"] >= 18,
                 "host_scans": host["successful_snapshots"] >= (3 if mode == "sampled" else 2),
                 "no_host_errors": host["observation_errors"] == 0,
                 "low_observed_competition": host["maximum_observed_foreign_average_cores"] < .25}
        if records != [source_record(path) for path in paths]:
            raise ValueError("Evidence changed during audit")
        rows.append({"index": index, "mode": mode, "scheduler": jobs[f"21647_{index}"],
                     "evidence": records, "replay": replay, "gates": gates,
                     "wall_s": measured["command_wall_s"], "cpu_s": summary["cpu_delta_usec"]["usage_usec"] / 1e6,
                     "mean_cpu_cores": summary["mean_cpu_cores_over_observed_span"],
                     "before_check_wall_s": verified["before_check_wall_s"],
                     "after_check_wall_s": verified["after_check_wall_s"],
                     "max_foreign_cores": host["maximum_observed_foreign_average_cores"],
                     "gnu_time": timing})
    paired = paired_summary(rows)
    return {**paired, "status": "verified_overhead_panel_evaluated", "runs": rows, "native_work": native_work,
            "all_protocol_gates_met": paired["numerical_wall_budget_met"] and all(all(r["gates"].values()) for r in rows),
            "before_after_identity_verified": True, "scientific_execution_authorized": False,
            "accounting": accounting, "protocol": protocol_record, "source": source_record(__file__),
            "limitations": ["Incremental periodic sampling cost, not comparison with an uninstrumented run.",
                            "A compute-only engineering fixture, not general overhead or scientific performance.",
                            "Full-file checks warm caches; their wall time is outside the command timer.",
                            "No overhead subtraction, selective exclusions or scientific launch authorization.",
                            "Before/after snapshots cannot prove absence of temporary changes or missed contention."]}


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

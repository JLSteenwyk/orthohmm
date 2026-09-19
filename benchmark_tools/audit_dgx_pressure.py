"""Describe retained host PSI counters; do not admit or correct native timing."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.summarize_frontier_overhead import terminal_scheduler_rows


def parse_pressure(raw, resource, *, system_scope=True):
    rows = {}
    for line in raw.splitlines():
        parts = line.split()
        if len(parts) != 5 or parts[0] not in ("some", "full") or parts[0] in rows:
            raise ValueError("Malformed or duplicate PSI row")
        values = {}
        for field in parts[1:]:
            key, separator, value = field.partition("=")
            if not separator or key in values:
                raise ValueError("Malformed or duplicate PSI field")
            values[key] = value
        if set(values) != {"avg10", "avg60", "avg300", "total"}:
            raise ValueError("Unexpected PSI fields")
        if not re.fullmatch(r"[0-9]+", values["total"]):
            raise ValueError("Invalid PSI total")
        if any(not math.isfinite(float(values[k])) or not 0 <= float(values[k]) <= 100
               for k in ("avg10", "avg60", "avg300")):
            raise ValueError("Invalid PSI average")
        rows[parts[0]] = int(values["total"])
    required = {"some"} if resource == "cpu" and system_scope else {"some", "full"}
    if resource not in ("cpu", "memory", "io") or not required <= rows.keys():
        raise ValueError("Missing PSI category or invalid resource")
    return {key: rows[key] for key in sorted(required)}


def pressure_summary(samples, resource):
    if len(samples) < 2:
        raise ValueError("Need at least two host observations")
    totals, identity, previous_end = [], None, -1
    field = f"host_{resource}_pressure"
    for sample in samples:
        start, end = sample["started_monotonic_ns"], sample["finished_monotonic_ns"]
        if any(type(v) is not int for v in (start, end)) or not previous_end < start <= end:
            raise ValueError("Invalid or overlapping observation windows")
        previous_end = end
        current = tuple(sample["raw"][key] for key in ("boot_id", "online_cpus", "cgroup_membership"))
        if identity is not None and identity != current:
            raise ValueError("Host or observer identity changed")
        identity = current
        if field not in sample["optional"] or any(e["field"] == field for e in sample["errors"]):
            raise ValueError("Missing or failed pressure read")
        value = parse_pressure(sample["optional"][field], resource)
        if totals and any(value[key] < totals[-1][key] for key in value):
            raise ValueError("Pressure counter decreased")
        totals.append(value)
    first, last = samples[0], samples[-1]
    span = (last["started_monotonic_ns"] + last["finished_monotonic_ns"]
            - first["started_monotonic_ns"] - first["finished_monotonic_ns"]) / 2e9
    delta = {key: totals[-1][key] - totals[0][key] for key in totals[0]}
    return dict(status="observed_pressure_totals", observations=len(samples), midpoint_span_s=span,
                stall_usec=delta, midpoint_percent={key: value / (span * 1e4) for key, value in delta.items()},
                cpu_full_undefined=(resource == "cpu"))


def read_record(path, records):
    raw = path.read_bytes()
    records.append(dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    return json.loads(raw)


def audit(archive):
    records = []
    accounting = archive / "accounting_terminal.txt"
    raw = accounting.read_bytes()
    scheduler = terminal_scheduler_rows(raw.decode(), 21838)
    records.append(dict(path=str(accounting), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    runs = []
    for index in range(18):
        directory = archive / "frontier_overhead_v1" / f"run_{index:02d}" / "measurement"
        samples = []
        points = sorted(directory.glob("point_*.json"))
        for path in points:
            point = read_record(path, records)
            if len(point["host"]) != 2:
                raise ValueError("Require paired host brackets")
            samples.extend(point["host"])
        row = dict(index=index, scheduler=scheduler[index], points=len(points), resources={})
        done_path = directory / "done.json"
        if done_path.exists():
            done = read_record(done_path, records)
            row["native_exit_code"] = done["exit_code"]
            row["recorded_window_encloses_native"] = bool(samples and
                samples[0]["finished_monotonic_ns"] < done["started_ns"]
                < done["finished_ns"] < samples[-1]["started_monotonic_ns"])
        else:
            row["recorded_window_encloses_native"] = False
        for resource in ("cpu", "memory", "io"):
            try:
                row["resources"][resource] = pressure_summary(samples, resource)
            except (ValueError, KeyError, TypeError) as error:
                row["resources"][resource] = dict(status="unavailable", reason=str(error))
        runs.append(row)
    for record in records:
        if hashlib.sha256(Path(record["path"]).read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError("Evidence changed during pressure audit")
    return dict(status="retained_host_pressure_diagnostic", runs=runs, records=records,
                scientific_timings_admitted=False, controlled_workload_verified=False,
                documentation="https://docs.kernel.org/accounting/psi.html",
                limitations=["Retrospective diagnostic with no new exclusion threshold or panel admission.",
                    "Host pressure includes native work and other tasks; it does not identify foreign interference.",
                    "CPU full is undefined at system scope and is not interpreted.",
                    "Midpoint percentages are descriptive; reads are non-atomic and pressure accounting may be delayed.",
                    "Zero recorded stalls do not prove isolation, no cache/I/O interference, or absence of throttling.",
                    "Incomplete native coverage and failed scheduler tasks remain visible; no missing values are zero-filled."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.archive.resolve())
    result["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")

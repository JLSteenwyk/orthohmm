"""Describe audited dual CPU flags; no causal attribution or eligibility changes."""

import argparse
import ast
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess

from benchmark_tools.measure_native_hierarchy_step import interval_point
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.summarize_pressure_panel_flags import distribution

AUDIT_SHA = "f6451d1fd5a96ff5f1a9b2155c6e4d9f17fabeb4de8a86cd53fbd1e70f37bff3"
CORE = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
STAGES = ("search", "edge_thresholds", "network_edges", "clustering", "profile_expansion",
          "refinement", "phylogeny_candidates", "phylogeny", "orthogroup_materialization")


def counter(snapshot, name):
    matches = [line.split() for line in snapshot["raw"]["proc_stat"].splitlines()
               if line.split() and line.split()[0] == name]
    if len(matches) != 1 or len(matches[0]) != 2 or not matches[0][1].isascii() or not matches[0][1].isdigit():
        raise ValueError("Missing or invalid proc counter: " + name)
    return int(matches[0][1])


def phase_bounds(metrics, native_wall):
    stages = metrics["stages"]
    if not stages or set(stages) - set(STAGES):
        raise ValueError("Unknown native metrics stage")
    if type(native_wall) not in (int, float) or not math.isfinite(native_wall) or native_wall <= 0:
        raise ValueError("Invalid enclosing command duration")
    durations = {key: value["wall_s"] for key, value in stages.items()}
    if any(type(v) not in (int, float) or not math.isfinite(v) or v < 0 for v in durations.values()):
        raise ValueError("Invalid stage duration")
    # Round-to-microsecond stage durations do not locate the uninstrumented gaps.
    error = len(stages) * 1e-6
    slack = native_wall - sum(durations.values())
    if slack < -error:
        raise ValueError("Stages exceed the enclosing native command")
    prefix, bounds = 0., []
    for name in STAGES:
        if name not in durations:
            continue
        duration = durations[name]
        bounds.append(dict(stage=name, earliest_start_s=max(0., prefix-error),
            latest_start_s=prefix+max(0., slack)+error,
            earliest_end_s=max(0., prefix+duration-error),
            latest_end_s=prefix+duration+max(0., slack)+error))
        prefix += duration
    return bounds


def certain_phase(bounds, start, end):
    phases = [b["stage"] for b in bounds if start >= b["latest_start_s"] and end <= b["earliest_end_s"]]
    if len(phases) > 1:
        raise ValueError("Overlapping certain stages")
    return phases[0] if phases else None


def summarize_rows(rows):
    fields = ("residual_cores", "native_average_cores", "overhang_s", "outside_frontier_cpu_s",
              "root_minus_frontier_cpu_s", "batch_step_cpu_s", "process_creations", "context_switches",
              "native_cpu_pressure_some_usec", "root_direct_processes")
    return dict(count=len(rows), distributions={k: distribution([r[k] for r in rows]) for k in fields} if rows else {},
                certain_phases=dict(Counter(r["certain_phase"] or "unresolved" for r in rows)))


def describe(run, measured, bounds):
    screen = run["screening"]
    if screen != measured["screening"]:
        raise ValueError("Audited and retained screening differ")
    points, narrow = measured["points"], screen["narrow_intervals"]
    original = screen["original_screening"]
    if (len(points) != len(narrow)+1 or any(len(original[k]) != len(narrow) for k in
            ("frontier_intervals", "hierarchy_intervals", "native_pressure_intervals"))):
        raise ValueError("Incomplete interval coverage")
    flags = [i for i, r in enumerate(narrow) if not r["screen_passed"]]
    if flags != run["narrow_flagged_intervals"] or flags != screen["narrow_flagged_intervals"]:
        raise ValueError("Audited narrow flags differ")
    rows = []
    for i, (left, right, value) in enumerate(zip(points, points[1:], narrow)):
        before, after = left["host"][0], right["hierarchy_host_after"]
        deltas = {key: counter(after, key)-counter(before, key) for key in ("processes", "ctxt")}
        if min(deltas.values()) < 0:
            raise ValueError("Host proc counter decreased")
        a, b = (interval_point(p, run["job_id"])["native_read_ns"] for p in (left, right))
        start, end = ((sum(t)/2-run["started_ns"])/1e9 for t in (a, b))
        if not math.isclose(end-start, value["wall_s"], abs_tol=1e-8):
            raise ValueError("Native interval timestamps disagree")
        frontier = original["frontier_intervals"][i]
        rows.append(dict(index=i, flagged=i in flags, reasons=value["reasons"],
            start_s=start, end_s=end, certain_phase=certain_phase(bounds, start, end),
            residual_cores=value["signed_unassigned_average_cores"],
            native_average_cores=value["native_cpu_s"]/value["wall_s"], overhang_s=value["outer_read_overhang_s"],
            outside_frontier_cpu_s=frontier["outside_target_frontier_cpu_s"],
            root_minus_frontier_cpu_s=frontier["root_minus_frontier_cpu_s"],
            batch_step_cpu_s=original["hierarchy_intervals"][i]["step_cpu_s"]["step_batch"],
            process_creations=deltas["processes"], context_switches=deltas["ctxt"],
            native_cpu_pressure_some_usec=original["native_pressure_intervals"][i]["native_stall_usec"]["cpu"]["some"],
            root_direct_processes=max(p["frontier"][k]["ancestor_direct_process_counts"]["/"]
                for p in (left, right) for k in ("inventory_before", "inventory_after")),
            active_outside_scopes={k:v for k,v in frontier["scope_cpu_s"].items()
                if k != left["frontier"]["target"] and v > 0}))
    phases = {}
    for name in sorted({r["certain_phase"] or "unresolved" for r in rows}):
        selected = [r for r in rows if (r["certain_phase"] or "unresolved") == name]
        phases[name] = dict(all=summarize_rows(selected),
            flagged=summarize_rows([r for r in selected if r["flagged"]]),
            unflagged=summarize_rows([r for r in selected if not r["flagged"]]))
    return dict(index=run["index"], method=run["method"], phase_bounds=bounds, intervals=rows,
                by_certain_phase=phases,
                all=summarize_rows(rows), flagged=summarize_rows([r for r in rows if r["flagged"]]),
                unflagged=summarize_rows([r for r in rows if not r["flagged"]]))


def report(audit_path, repo):
    source = record(audit_path)
    if source["sha256"] != AUDIT_SHA:
        raise ValueError("Wrong complete audit")
    audit = json.loads(gzip.decompress(audit_path.read_bytes()))
    if [r["index"] for r in audit["runs"]] != [0, 1, 2] or audit["validated_tasks"] != 3:
        raise ValueError("Require complete validated native panel")
    code = subprocess.check_output(["git", "show", CORE+":orthohmm/orthohmm.py"], cwd=repo)
    calls = sorted((n for n in ast.walk(ast.parse(code)) if isinstance(n, ast.Call)
        and isinstance(n.func, ast.Attribute) and isinstance(n.func.value, ast.Name)
        and n.func.value.id == "metrics" and n.func.attr == "stage"), key=lambda n:n.lineno)
    order = list(dict.fromkeys(n.args[0].value for n in calls))
    if order != [*STAGES, "output"]:
        raise ValueError("Frozen source stage order changed")
    rows, evidence = [], [source]
    for run in audit["runs"]:
        if run["status"] != "validated":
            raise ValueError("Cannot drop an invalid task")
        reports = [r for r in run["inventory"] if Path(r["path"]).name == "dual_bracket_report.json"]
        if len(reports) != 1:
            raise ValueError("Ambiguous retained report")
        item = reports[0]
        check(item)
        measured = json.loads(Path(item["path"]).read_text())
        evidence.append(item)
        bounds = []
        if run["method"].startswith("orthohmm_"):
            items = [r for r in run["inventory"] if Path(r["path"]).name == run["method"]+".json"]
            if len(items) != 1:
                raise ValueError("Missing native metrics")
            check(items[0])
            evidence.extend(items)
            bounds = phase_bounds(json.loads(Path(items[0]["path"]).read_text()), run["native_wall_s"])
        rows.append(describe(run, measured, bounds))
    for item in evidence:
        check(item)
    return dict(status="dual_native_flags_described_not_causally_attributed", runs=rows, evidence=evidence,
        frozen_stage_source=dict(commit=CORE, path="orthohmm/orthohmm.py", bytes=len(code), sha256=hashlib.sha256(code).hexdigest()),
        source=record(__file__), scientific_timings_admitted=False, publication_ready=False,
        limitations=["Post-outcome exploratory description, not new raw replay, acceptance rules or causal inference.",
            "Host, hierarchy, pressure and frontier windows differ; their values cannot be causally subtracted.",
            "Host process creation includes all tasks, not just native subprocesses.",
            "Stage bounds allocate all uninstrumented time conservatively; boundary intervals remain unresolved.",
            "OrthoFinder phase assignment is unavailable in this analysis.",
            "Root-direct task counts do not identify which tasks used CPU; transient tasks and accounting delay remain unresolved.",
            "No flags removed, thresholds changed, selective reruns or scientific timing admission."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = report(args.audit.resolve(), args.repo.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")

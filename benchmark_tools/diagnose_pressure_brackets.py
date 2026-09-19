"""Compare retained host brackets without changing historical timing eligibility."""

import argparse
from collections import Counter
import json
from pathlib import Path

from benchmark_tools.measure_native_hierarchy_step import interval_point
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_host_counters import summarize as host_delta
from benchmark_tools.probe_interval_cpu import interval
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.summarize_pressure_panel_flags import distribution


def compare(report):
    points, job = report["points"], report["job_id"]
    if len(points) < 2:
        raise ValueError("Require at least two observation points")
    outer, narrow = [], []
    for point in points:
        earlier, later = point["hierarchy_host_after"], point["host"][1]
        host_delta(earlier, later, point["ticks"])
        outer.append(interval_point(point, job))
        narrow.append(interval_point({**point, "host": [point["host"][0], earlier]}, job))
    broad_rows = [interval(a, b, job) for a, b in zip(outer, outer[1:])]
    original = report["screening"]["original_threshold_screen"]
    if broad_rows != original["intervals"]:
        raise ValueError("Outer replay differs from audited interval results")
    flags = [i for i, row in enumerate(broad_rows) if not row["screen_passed"]]
    if flags != original["flagged_intervals"]:
        raise ValueError("Outer flag inventory differs")
    narrow_rows = [interval(a, b, job) for a, b in zip(narrow, narrow[1:])]
    transitions, extra_host, extra_span = Counter(), [], []
    for broad, short in zip(broad_rows, narrow_rows):
        if broad["native_cpu_s"] != short["native_cpu_s"] or broad["wall_s"] != short["wall_s"]:
            raise ValueError("Native counter or duration changed")
        transitions[f"{broad['screen_passed']}->{short['screen_passed']}"] += 1
        extra_host.append(broad["host_busy_cpu_s"] - short["host_busy_cpu_s"])
        extra_span.append(broad["outer_read_overhang_s"] - short["outer_read_overhang_s"])
    return dict(intervals=len(broad_rows), transitions=dict(transitions),
                outer_flagged=len(flags),
                narrow_flagged=sum(not r["screen_passed"] for r in narrow_rows),
                narrow_reasons=dict(Counter(reason for r in narrow_rows for reason in r["reasons"])),
                added_host_cpu_s=distribution(extra_host), added_read_overhang_s=distribution(extra_span),
                outer_residual_cores=distribution([r["signed_unassigned_average_cores"] for r in broad_rows]),
                narrow_residual_cores=distribution([r["signed_unassigned_average_cores"] for r in narrow_rows]),
                narrow_overhang_s=distribution([r["outer_read_overhang_s"] for r in narrow_rows]),
                scientific_timings_admitted=False)


def run(source_path, digest, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    source = record(source_path)
    summary = read_frozen(source_path, digest)
    if (summary["status"] != "audited_pressure_flags_described"
            or [r["index"] for r in summary["runs"]] != list(range(18))):
        raise ValueError("Require complete audited panel summary")
    reports = {str(Path(item["path"]).parent.parent): item for item in summary["evidence"]
               if Path(item["path"]).name == "frontier_report.json"}
    rows, evidence = [], [source]
    for row in summary["runs"]:
        result = {k: row[k] for k in ("index", "method", "mode", "status")}
        if row["status"] == "validated" and row["mode"] == "periodic":
            matches = [item for directory, item in reports.items()
                       if Path(directory).name == f"run_{row['index']:02d}"]
            if len(matches) != 1:
                raise ValueError("Require exactly one audited task report")
            item = matches[0]
            check(item)
            result["bracket_diagnostic"] = compare(json.loads(Path(item["path"]).read_text()))
            if result["bracket_diagnostic"]["outer_flagged"] != row["diagnostic"]["flagged_intervals"]:
                raise ValueError("Flag summary differs from raw replay")
            evidence.append(item)
        rows.append(result)
    for item in evidence:
        check(item)
    result = dict(status="retained_pressure_brackets_compared", runs=rows, evidence=evidence,
                  source=record(__file__), scientific_timings_admitted=False,
                  limitations=["Retrospective diagnostic, not a new eligibility rule or complete frontier re-audit.",
                               "Both bracket choices retain non-atomic counter reads and accounting delay.",
                               "Identical native counters isolate the arithmetic effect of the later host read, not causal runtime interference.",
                               "Failed tasks remain failed; no wall-time correction or scientific timing admission."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.summary.resolve(), args.sha256, args.output.absolute())

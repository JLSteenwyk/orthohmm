"""Endpoint-recomputed accounting at fixed scales; never timing eligibility."""

import argparse
import gzip
import json
from pathlib import Path

from benchmark_tools.probe_dual_cpu_brackets import compare
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.summarize_pressure_panel_flags import distribution

NATIVE_SHA = "f6451d1fd5a96ff5f1a9b2155c6e4d9f17fabeb4de8a86cd53fbd1e70f37bff3"
CONTROL_SHA = "bd8e11ead6ee2e979784c8cd97730d3cf82ee62e4343c9cdf18c0475d2221cd3"
SCALES = (5, 10, 30)


def windows(point_count, width):
    if type(point_count) is not int or point_count < 2 or type(width) is not int or width <= 0:
        raise ValueError("Invalid endpoint window inventory")
    return [(start, min(start+width, point_count-1)) for start in range(0, point_count-1, width)]


def describe(measured):
    points, job = measured["points"], measured["job_id"]
    flags = measured["screening"]["narrow_flagged_intervals"]
    if len(points) != len(measured["screening"]["narrow_intervals"])+1:
        raise ValueError("Incomplete one-sample interval inventory")
    if flags != [i for i, row in enumerate(measured["screening"]["narrow_intervals"]) if not row["screen_passed"]]:
        raise ValueError("Original flags disagree")
    scales = []
    for width in SCALES:
        rows = []
        for start, end in windows(len(points), width):
            result = compare(points[start], points[end], job, enforce_gap=False)
            rows.append(dict(start_index=start, end_index=end, partial=end-start < width,
                original_flag_indices=[i for i in flags if start <= i < end], brackets=result))
        scales.append(dict(width_in_sample_steps=width, windows=rows,
            partial_windows=sum(row["partial"] for row in rows),
            narrow_flagged_windows=sum(not row["brackets"]["narrow"]["screen_passed"] for row in rows),
            residual_cores=distribution([row["brackets"]["narrow"]["signed_unassigned_average_cores"] for row in rows])))
    whole = compare(points[0], points[-1], job, enforce_gap=False)
    if whole != measured["screening"]["observation_window"]:
        raise ValueError("Whole observation window does not reproduce")
    return dict(original_interval_count=len(points)-1, original_flag_indices=flags,
                scales=scales, whole_observation_window=whole,
                whole_command_original_screen=measured["screening"]["original_screening"]["original_threshold_screen"]["whole_command_screen"])


def load(path, sha):
    source = record(path)
    if source["sha256"] != sha:
        raise ValueError("Wrong audited source")
    result = json.loads(gzip.decompress(path.read_bytes()))
    check(source)
    return result, source


def report(native_path, control_path):
    native, native_source = load(native_path, NATIVE_SHA)
    controls, control_source = load(control_path, CONTROL_SHA)
    if native["validated_tasks"] != 3 or controls["validated_trials"] != 9:
        raise ValueError("Require both complete validated panels")
    rows, evidence = [], [native_source, control_source]
    for run in native["runs"]:
        if run["status"] != "validated":
            raise ValueError("Cannot omit failed native diagnostic")
        files = [item for item in run["inventory"] if Path(item["path"]).name == "dual_bracket_report.json"]
        if len(files) != 1:
            raise ValueError("Missing native measurement")
        check(files[0])
        measured = json.loads(Path(files[0]["path"]).read_text())
        evidence.extend(files)
        rows.append(dict(kind="native", index=run["index"], method=run["method"], accounting=describe(measured)))
    for run in controls["trials"]:
        if run["status"] != "validated":
            raise ValueError("Cannot omit failed control")
        rows.append(dict(kind="control", index=run["index"], method=run["mode"],
                         accounting=describe(run["replay"]["measurement"]["measured"])))
    for item in evidence:
        check(item)
    return dict(status="cpu_window_scales_described", runs=rows, evidence=evidence, source=record(__file__),
        scientific_timings_admitted=False, publication_ready=False,
        limitations=["Post-outcome accounting diagnostic at fixed sample-index widths, not prespecified scientific endpoints.",
        "Every final partial window is retained; widths in samples are not exactly elapsed seconds.",
        "Endpoint counters are recomputed, never sums of overlapping residual windows.",
        "Longer windows can conceal brief real interference as well as average non-atomic accounting noise.",
        "Existing absolute/relative screen checks are descriptive at these scales, not calibrated new acceptance rules.",
        "Whole observation windows include wrapper work and are not exact command windows.",
        "All original flags retained; no threshold changes, causal attribution or timing admission."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("native-audit", "control-audit", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    result = report(args.native_audit.resolve(), args.control_audit.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")

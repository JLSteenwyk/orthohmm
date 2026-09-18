"""Decompose retained CPU-counter windows without changing timing admission."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.measure_native_interval_step import evaluate
from benchmark_tools.probe_host_counters import summarize
from benchmark_tools.probe_interval_cpu import interval
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

INPUTS = (
    "6e76afc4379ca894fecbec61bc2fb1c44de408fce6e4c2feefe21a8e94f13a46",
    "2d63f27662ec8aadd7a313e0ba42617f52f4d8caa81a4a3e3c87fb67a56560e3",
    "30daef66f09b6f5c07c12fdac6caa3ce5d176cbd4a6e234accb3a5a74113099d",
)
METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full")


def parse_stat(raw):
    fields = {}
    for line in raw.splitlines():
        key, value = line.split()
        if key in fields:
            raise ValueError("Duplicate CPU statistic")
        fields[key] = int(value)
        if fields[key] < 0:
            raise ValueError("Negative CPU statistic")
    if not {"usage_usec", "user_usec", "system_usec"} <= fields.keys():
        raise ValueError("Missing CPU accounting fields")
    return fields


def delta_stat(left, right):
    a, b = parse_stat(left), parse_stat(right)
    result = {key: (b[key] - a[key]) / 1e6 for key in ("usage_usec", "user_usec", "system_usec")}
    if any(v < 0 for v in result.values()):
        raise ValueError("CPU accounting decreased")
    return result


def decompose(left, right, job):
    original = interval(left, right, job)
    ticks = left["ticks_per_second"]
    outer = summarize(left["host"][0], right["host"][1], ticks)
    inner = summarize(left["host"][1], right["host"][0], ticks)
    left_window = summarize(*left["host"], ticks)["accounted_host_busy_cpu_s"]
    right_window = summarize(*right["host"], ticks)["accounted_host_busy_cpu_s"]
    mass = outer["accounted_host_busy_cpu_s"] - inner["accounted_host_busy_cpu_s"]
    if abs(mass - left_window - right_window) > 1e-9:
        raise ValueError("Host window telescoping identity failed")
    native = delta_stat(left["native_cpu_stat"], right["native_cpu_stat"])
    observer = delta_stat(left["host"][0]["optional"]["cgroup_cpu.stat"],
                          right["host"][1]["optional"]["cgroup_cpu.stat"])
    return dict(original_screen=original, host_outer_busy_cpu_s=outer["accounted_host_busy_cpu_s"],
        host_inner_busy_cpu_s=inner["accounted_host_busy_cpu_s"], read_window_busy_cpu_s=mass,
        host_outer_field_cpu_s={k: v / ticks for k, v in outer["cpu_delta_ticks"].items()},
        native_step_delta_cpu_s=native, observer_leaf_delta_cpu_s=observer,
        inner_host_minus_native_cpu_s=inner["accounted_host_busy_cpu_s"] - native["usage_usec"],
        outer_host_minus_native_minus_observer_cpu_s=(original["signed_unassigned_cpu_s"] - observer["usage_usec"]),
        controlled_workload_verified=False, scientific_timings_admitted=False)


def diagnose(report):
    replay = evaluate(report["points"], report["native"], report["job_id"])
    if replay != report["screening"]:
        raise ValueError("Stored screening differs from raw counter replay")
    rows = [dict(index=i, **decompose(a, b, report["job_id"]))
            for i, (a, b) in enumerate(zip(report["points"], report["points"][1:]))]
    return dict(job_id=report["job_id"], original_flagged_intervals=replay["flagged_intervals"],
                original_whole_command_screen=replay["whole_command_screen"], intervals=rows)


def run(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    base = root / "benchmarks/work/dgx_interval_native_21810/interval_native_smoke_v2"
    records, rows = [], []
    for i, (method, sha) in enumerate(zip(METHODS, INPUTS)):
        path = base / f"run_{i:02d}/measurement/interval_report.json"
        identity = record(path)
        if identity["sha256"] != sha:
            raise ValueError("Changed retained native interval evidence")
        rows.append(dict(method=method, evidence=identity, **diagnose(json.loads(path.read_text()))))
        records.append(identity)
    for item in records:
        check(item)
    result = dict(status="retrospective_native_interval_counter_decomposition", runs=rows,
        source=record(__file__), helpers=[record(Path(__file__).with_name(n)) for n in (
            "measure_native_interval_step.py", "probe_host_counters.py", "probe_interval_cpu.py",
            "screen_bracketed_cpu.py", "probe_dgx_step_separation.py")],
        publication_ready=False, controlled_workload_verified=False, scientific_timings_admitted=False,
        limitations=["Retrospective diagnostic on all three retained smokes; no threshold or inclusion rule changes.",
            "Inner/outer counter differences describe reported accounting, not actual CPU-time bounds or accounting delay.",
            "Observer leaf excludes other batch-step/system work; it is not a complete observer-overhead measurement.",
            "Residual arithmetic is not attribution to unrelated work and must not correct native wall time.",
            "Overlapping host windows must not be summed; negative residuals are preserved.",
            "User/system splits differ in accounting scope and granularity; non-CPU interference is unmeasured."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.absolute())

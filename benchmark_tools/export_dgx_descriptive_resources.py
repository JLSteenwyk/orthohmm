"""Export all DGX resource observations without admitting controlled comparisons."""

import argparse
import csv
import json
from pathlib import Path
import statistics
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_sequence_numeric import completed
from benchmark_tools.audit_dgx_scientific_metadata import PLAN_SHA
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def assemble(native, resource, host, plan):
    if (native["status"] != "native_panel_validation_complete_not_timing_admission"
            or resource["failures"] != 0 or host["review_failures"] != 0
            or any(len(item["runs"]) != 27 for item in (native, resource, host, plan))):
        raise ValueError("Require complete reviewed panel")
    rows = []
    for i, (n, r, h, p) in enumerate(zip(native["runs"], resource["runs"], host["runs"], plan["runs"])):
        if (any(item["index"] != i for item in (n, r, h, p))
                or any(n[k] != h[k] for k in ("method", "proteomes", "repeat"))
                or n["method"] != p["native_method"] or n["proteomes"] != p["proteomes"] or n["repeat"] != p["repeat"]
                or r["job_id"] != h["host_review"]["job_id"]
                or r["status"] != "resource_accounting_reproduced_not_timing_admission"
                or n["status"] not in ("native_outputs_validated", "native_validation_failed")):
            raise ValueError("Resource run identities differ")
        timing, summary = r["gnu_time"], r["summary"]
        workload = h["host_review"]["retained_host_summary"]
        rows.append({"index": i, "method": n["method"], "proteomes": n["proteomes"], "repeat": n["repeat"],
            "proteins": p["dataset"]["proteins"], "sequence_characters": p["dataset"]["sequence_characters"],
            "native_validation_status": n["status"], "native_validation_error": n.get("error", ""),
            "native_elapsed_seconds": timing["elapsed_seconds"],
            "native_cpu_seconds": timing["user_seconds"] + timing["system_seconds"],
            "gnu_max_process_rss_kib": timing["max_process_rss_kib"],
            "cgroup_peak_bytes": summary["maximum_reported_cgroup_peak_bytes"],
            "sampled_aggregate_rss_max_bytes": summary["maximum_sampled_sum_rss_bytes"],
            "resource_observations": summary["observations"],
            "resource_samples_with_process_errors": summary["samples_with_process_errors"],
            "host_status": workload["status"], "scientific_timing_admitted": False})
    if sum(r["native_validation_status"] == "native_validation_failed" for r in rows) != native["failures"]:
        raise ValueError("Native failure total differs")
    summaries = []
    for method, size in sorted({(r["method"], r["proteomes"]) for r in rows}):
        group = [r for r in rows if (r["method"], r["proteomes"]) == (method, size)]
        if len(group) != 3 or {r["repeat"] for r in group} != {0, 1, 2}:
            raise ValueError("Incomplete repeated-run cell")
        valid = all(r["native_validation_status"] == "native_outputs_validated" for r in group)
        result = {"method": method, "proteomes": size, "runs": 3,
            "native_valid_runs": sum(r["native_validation_status"] == "native_outputs_validated" for r in group),
            "host_inconclusive_runs": sum(r["host_status"] == "inconclusive" for r in group),
            "rss_incomplete_runs": sum(r["resource_samples_with_process_errors"] > 0 for r in group),
            "status": "descriptive_only" if valid else "native_validation_incomplete", "scientific_timing_admitted": False}
        for metric in ("native_elapsed_seconds", "native_cpu_seconds", "gnu_max_process_rss_kib", "cgroup_peak_bytes",
                       "sampled_aggregate_rss_max_bytes"):
            values = [r[metric] for r in group]
            result[metric] = {"median": statistics.median(values), "minimum": min(values), "maximum": max(values)} if valid else None
        summaries.append(result)
    return {"status": "descriptive_resource_observations_not_controlled_comparison", "runs": rows, "summaries": summaries,
        "scientific_timings_admitted": 0, "publication_ready": False,
        "limitations": ["All repeats retained; no fastest-run selection, speedup ratios or complexity fits.",
            "Cells with native validation failures retain observations but have null aggregate metrics.",
            "Host classifications remain unchanged; medians do not resolve isolation uncertainty.",
            "Sampled aggregate RSS is incomplete where process reads failed and can miss between-sample peaks.",
            "GNU-time maximum process RSS, sampled aggregate RSS and cgroup memory peak are not interchangeable.",
            "Cgroup peak includes charged cache/kernel memory and may predate native launch.",
            "Only three repetitions on one nested proteome series; taxa and dataset size co-vary."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--job", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    scheduler = completed(args.job, 2, "64G")
    root = args.root.resolve()
    results = root / "benchmark_tools/results"
    paths = [root / "benchmarks/work/dgx_native_output_validation_20260918.json",
             results / "dgx_resource_replay_20260918.json", results / "dgx_completed_panel_review_20260918.json"]
    inputs = [record(p) for p in paths]
    native, resource, host = [read_pinned(p, r["sha256"]) for p, r in zip(paths, inputs)]
    executor = root / "benchmarks/work/publication_dgx_native_validation_v1/benchmark_tools/validate_dgx_native_panel.py"
    if native["source"] != record(executor) or native["metadata_review"] != inputs[2]:
        raise ValueError("Native validator source or metadata binding differs")
    check(native["archive_inventory"])
    plan_path = results / "dgx_scaling_commands_20260917.json"
    report = assemble(native, resource, host, read_pinned(plan_path, PLAN_SHA))
    checked = [*inputs, native["archive_inventory"], native["source"], *native["helpers"], record(plan_path)]
    for row in native["runs"]:
        checked.extend(row.get("validation", {}).get("checked_files", []))
    for row in resource["runs"]:
        checked.extend(row["evidence"])
    for item in checked:
        check(item)
    args.output.mkdir(parents=True)
    report.update(source=record(__file__), inputs=inputs, checked_records=checked, scheduler=scheduler)
    with (args.output / "runs.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(report["runs"][0]))
        writer.writeheader()
        writer.writerows(report["runs"])
    with (args.output / "report.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")

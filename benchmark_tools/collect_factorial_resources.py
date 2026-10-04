"""Consolidate retained factorial stage costs without measuring or imputing new costs."""

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path


CELLS = tuple("p%d_c%d_r%d" % factors for factors in itertools.product(range(2), repeat=3))
ARMS = tuple("p%d_c%d" % factors for factors in itertools.product(range(2), repeat=2))
FIELDS = ("wall_s", "user_cpu_s", "system_cpu_s", "peak_process_tree_rss_bytes")
FIXED_INPUTS = {
    "ob_results": ("benchmark_tools/results/orthobench_factorial_results_20260916.json",
                   "6a0d588b5cb47c60fc6bc8bae8aa0c83e5f2aadb11de970919d8c6527c387141"),
    "ob_preparation": ("benchmark_tools/results/orthobench_factorial_prepared_20260916.json",
                       "5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382"),
    "qfo_preparation": ("benchmarks/results/qfo_corrected_factorial_v1/manifest.json",
                        "d8385c50426e690afd6d32f3c5302e678de6977c9841451d013201b0f75b564a"),
}
QFO_ADMISSIONS = {
    "p0_c0_r1": (21761, "97341e1b9ef6ac36e6b8c329aa3b5161690ccf00e617895f1d4ccb6ea127b16a"),
    "p0_c1_r1": (21762, "27d359b059d3c7589e8f85adf03a00fc5e53d57e6b80f4e8b75727f0749841a8"),
    "p1_c0_r1": (21763, "36f3c8a3b9cebf8f474b3cb956d3e3e79cebb5243a4baa52e89a1d0a5c23d317"),
    "p1_c1_r1": (21764, "c8a289ac8128711da6c4e93654ce6953e48854b8d21b2adfea074d3feb8c70aa"),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def record(root, path):
    path = Path(path)
    if not path.is_absolute():
        path = root / path
    path = path.resolve()
    require(path.is_relative_to(root), "Input outside repository")
    data = path.read_bytes()
    return {"path": str(path.relative_to(root)), "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest()}


def load(root, reference):
    observed = record(root, reference["path"])
    require(all(observed[key] == reference[key] for key in ("bytes", "sha256")),
            "Retained input checksum mismatch")
    return json.loads((root / observed["path"]).read_text()), observed


def finite(value, name, positive=False):
    require(not isinstance(value, bool) and isinstance(value, (int, float))
            and math.isfinite(value) and (value > 0 if positive else value >= 0),
            "Invalid resource value: " + name)
    return value


def validate_plan(plan):
    cells = plan["cells"]
    require(len(cells) == 8 and {c["label"] for c in cells} == set(CELLS),
            "Incomplete or duplicate factorial cells")
    require(set(plan["candidate_arms"]) == set(ARMS), "Incomplete candidate arms")
    for cell in cells:
        label = cell["label"]
        expected = (label[1] == "1", label[4] == "1", label[7] == "1")
        observed = tuple(cell[k] for k in ("profile_expansion", "candidate_expansion", "reconciliation"))
        require(all(type(v) is bool for v in observed) and observed == expected,
                "Factor flags differ from cell identity")
        require(cell["runtime_kind"] == "incremental_cached_replay",
                "Unexpected factorial runtime scope")
    for arm in plan["candidate_arms"].values():
        finite(arm["incremental_preparation_seconds"], "candidate preparation", positive=True)


def metrics_values(metrics):
    require(metrics["status"] == "complete", "Incomplete native metrics")
    require(metrics["metadata"]["cpu_budget"] == 32, "Unexpected recorded CPU budget")
    require(metrics["rss_measurement"] == "sampled_sum_of_linux_proc_tree_rss",
            "Unexpected historical memory convention")
    values = {key: finite(metrics[key], key, positive=key in ("wall_s", "peak_process_tree_rss_bytes"))
              for key in FIELDS}
    require(type(values["peak_process_tree_rss_bytes"]) is int, "RSS bytes must be an integer")
    return values


def project(dataset, plan, measurements, batches):
    validate_plan(plan)
    require(set(measurements) == {c for c in CELLS if c.endswith("r1")},
            "Incomplete or unexpected reconciliation resources")
    require(set(batches) == set(measurements), "Incomplete original batch identities")
    rows = []
    for label in CELLS:
        arm = label[:-3]
        costs = metrics_values(measurements[label]) if label in measurements else dict.fromkeys(FIELDS)
        rows.append({
            "dataset": dataset, "cell": label, "candidate_preparation_arm": arm,
            "candidate_preparation_wall_s": plan["candidate_arms"][arm]["incremental_preparation_seconds"],
            "candidate_preparation_cpu_s": None, "candidate_preparation_peak_bytes": None,
            "candidate_preparation_scope": "Partition copy/validation and optional satellite expansion from cached seeds; shared between R cells",
            "reconciliation_wall_s": costs["wall_s"],
            "reconciliation_user_cpu_s": costs["user_cpu_s"],
            "reconciliation_system_cpu_s": costs["system_cpu_s"],
            "reconciliation_peak_sampled_tree_rss_bytes": costs["peak_process_tree_rss_bytes"],
            "reconciliation_mean_cpu_cores": None if costs["wall_s"] is None else (
                costs["user_cpu_s"] + costs["system_cpu_s"]) / costs["wall_s"],
            "reconciliation_measurement_status": "recorded_incremental" if label in measurements else "not_applicable",
            "reconciliation_scope": "Cached-candidate native reconciliation; excludes initial search, profile/seed construction, preparation, conversion and scoring",
            "original_batch": batches.get(label),
            "full_pipeline_wall_s": None, "full_pipeline_cpu_s": None,
            "full_pipeline_peak_memory_bytes": None,
            "full_pipeline_measurement_status": "unavailable_not_measured_per_cell",
            "unknown_shared_host_contention": True,
        })
    return rows


def collect(root):
    root = Path(root).resolve()
    inputs, documents = [], {}
    for key, (name, digest) in FIXED_INPUTS.items():
        ref = record(root, name)
        require(ref["sha256"] == digest, "Unexpected fixed factorial input: " + key)
        documents[key], ref = load(root, ref)
        inputs.append(ref)
    ob = documents["ob_results"]
    measurements, batches = {}, {}
    for label, admission in ob["native_validation"].items():
        require(admission["status"] == "native_group_output_verified", "OrthoBench cell not admitted")
        require(admission["cell"] == label, "OrthoBench cell identity differs")
        measurements[label], ref = load(root, admission["native_metrics"])
        inputs.append(ref)
        batches[label] = admission["integrity"]["scheduler"]
        values = metrics_values(measurements[label])
        require(all(ob["coverage_resources"][label][key] == values[key] for key in FIELDS),
                "OrthoBench resource table differs from native metrics")
    rows = project("OrthoBench", documents["ob_preparation"], measurements, batches)
    measurements, batches = {}, {}
    for label, (job, digest) in QFO_ADMISSIONS.items():
        name = "benchmark_tools/results/qfo_corrected_factorial_native_admission_%d.json" % job
        ref = record(root, name)
        require(ref["sha256"] == digest, "Unexpected corrected QfO admission")
        admission, ref = load(root, ref)
        inputs.append(ref)
        require(admission["status"] == "corrected_qfo_native_pair_output_verified"
                and admission["cell"] == label, "Wrong corrected QfO cell")
        native = admission["native_group_integrity"]
        require(native["status"] == "native_group_output_verified" and native["cell"] == label,
                "Wrong corrected QfO native group admission")
        measurements[label], ref = load(root, native["native_metrics"])
        inputs.append(ref)
        metrics = measurements[label]
        arm = documents["qfo_preparation"]["candidate_arms"][label[:-3]]
        require(metrics["input"]["candidate_clusters"] == arm["candidate_partition"],
                "QfO metrics use different candidates")
        require(metrics["input"]["membership_constraints"] == arm.get("membership_constraints"),
                "QfO metrics use different constraints")
        batches[label] = native["integrity"]["scheduler"]
    rows.extend(project("Corrected QfO", documents["qfo_preparation"], measurements, batches))
    for ref in inputs:
        require(record(root, ref["path"]) == ref, "Input changed during collection")
    return {
        "schema": "factorial_retained_resources_v1", "status": "retained_stage_costs_consolidated",
        "inputs": inputs, "rows": rows, "cells": 16, "shared_candidate_arms": 8,
        "recorded_reconciliation_cells": 8, "full_pipeline_cost_cells_available": 0,
        "native_inference_or_scoring_repeated": False, "publication_ready": False,
        "limitations": [
            "Historical shared-host stage observations, not isolated tool speed or causal component overhead.",
            "Preparation intervals omit cache loading and seed/profile construction. An arm is shared by its R-off/R-on rows; do not sum them as independent preparation runs.",
            "Eight R-off reconciliation stages are not applicable, not zero-cost full pipelines.",
            "Peak is sampled sum of process-tree RSS, potentially double-counting shared pages and missing between-sample or short-lived peaks. It is not the newer scaling panel's lifetime cgroup peak.",
            "No observed contention series or newer collector-validity pass is imputed to historical runs.",
            "All sixteen per-cell full-pipeline costs are unavailable. Stage times are not added or substituted for them; stage peaks are neither added nor promoted to whole-run peaks.",
            "Original OrthoBench batch failures remain failed; separate scientific recovery does not turn those batches into clean completed runs.",
            "Direct source/admission/metrics pins are checked, not every transitive raw artifact or historical process-accounting assumption.",
        ],
    }


def render(report):
    lines = ["# Retained Factorial Stage Costs", "",
             "P: profile expansion; C: candidate expansion; R: reconciliation. P-off retains initial HMM search.",
             "All values are historical shared-host observations, not isolated efficiency comparisons.",
             "Preparation repeats the same shared arm in two rows; it excludes search/profile/seed construction.",
             "NA is unavailable/not applicable, not zero. RSS is sampled summed process-tree RSS, not lifetime peak.",
             "", "| Dataset | Cell | Shared preparation (s) | Reconciliation (s) | Mean CPU cores | Sampled tree RSS (GiB) | Full pipeline |",
             "| --- | --- | ---: | ---: | ---: | ---: | --- |"]
    for row in report["rows"]:
        def fmt(key, divisor=1, places=3):
            value = row[key]
            return "NA" if value is None else ("%." + str(places) + "f") % (value / divisor)
        lines.append("| %s | %s | %s | %s | %s | %s | NA |" % (
            row["dataset"], row["cell"], fmt("candidate_preparation_wall_s"),
            fmt("reconciliation_wall_s"), fmt("reconciliation_mean_cpu_cores"),
            fmt("reconciliation_peak_sampled_tree_rss_bytes", 2**30)))
    lines += ["", "## Limits", "", *["- " + text for text in report["limitations"]]]
    return "\n".join(lines) + "\n"


def write(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError("Refusing existing output directory")
    report = collect(root)
    report["source"] = record(root, Path(__file__))
    output.mkdir(parents=True)
    with (output / "resources.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    with (output / "resources.md").open("x") as stream:
        stream.write(render(report))
    columns = list(report["rows"][0])
    with (output / "resources.tsv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows({k: ("NA" if v is None else json.dumps(v, sort_keys=True)
                             if isinstance(v, dict) else v) for k, v in row.items()}
                        for row in report["rows"])
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    report = write(args.root, args.output_directory)
    print(json.dumps({k: report[k] for k in ("status", "cells", "recorded_reconciliation_cells",
                                            "full_pipeline_cost_cells_available")}))


if __name__ == "__main__":
    main()

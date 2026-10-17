"""Export retained resource intervals without aggregating incompatible scopes."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.summarize_matched_resources import verbose_time

DATASETS = ("OrthoBench", "QfO", "ThreeKingdoms")
EXCLUDE = {"command", "log", "time", "source", "inputs", "input", "input_records", "outputs",
           "runtime", "runtime_records", "metadata", "parameters", "counts", "semantics", "evidence",
           "retained_measurements", "scheduler", "memory", "host_before", "host_after"}
SCALAR_WALL = {"total_incremental_wall_s", "total_full_pipeline_and_scoring_wall_s", "total_pilot_wall_s",
               "native_and_scoring_wall_s", "candidate_wall_s", "scoring_recovery_wall_s"}
COLUMNS = ("method", "label", "dataset", "cell_index", "cell_status", "origin", "source_document",
           "source_location", "wall_seconds", "user_cpu_seconds", "system_cpu_seconds",
           "total_cpu_seconds_reported", "peak_memory_value", "peak_memory_unit", "memory_scope",
           "requested_cpu_budget", "allocated_cpus", "job_id", "measurement_scope", "cell_scope")
COLUMNS += ("full_inference", "cpu_scope", "measurement_status", "command_returncode", "gnu_exit_status",
            "scheduler_state", "scheduler_exit_code", "reference_location")


def number(value):
    if value is not None and (type(value) not in (int, float) or not math.isfinite(value) or value < 0):
        raise ValueError("Resource measurements must be finite nonnegative numbers or missing")
    return value


def first(node, names):
    return next((number(node[name]) for name in names if name in node), None)


def memory(node, context):
    explicit = node.get("memory")
    if isinstance(explicit, dict):
        return number(explicit.get("value")), explicit.get("unit"), explicit.get("scope")
    for key, unit, scope in (
        ("max_process_rss_kib", "KiB", "Reported maximum process RSS; not summed process-tree memory"),
        ("max_reported_process_rss_kib", "KiB", "Reported maximum process RSS; not summed process-tree memory"),
        ("peak_process_tree_rss_bytes", "bytes", context.get("rss_measurement", "Sampled process-tree RSS")),
        ("peak_process_tree_rss_gib", "GiB", context.get("rss_measurement", "Sampled process-tree RSS")),
        ("peak_process_rss_GiB", "GiB", "Reported process RSS"),
        ("parent_peak_rss_GiB", "GiB", "Reported parent-process RSS only"),
    ):
        if key in node:
            return number(node[key]), unit, node.get("memory_scope", scope)
    return None, None, None


def intervals(node, location=(), context=None):
    context = dict(context or {})
    if isinstance(node, list):
        for index, child in enumerate(node):
            yield from intervals(child, (*location, index), context)
        return
    if not isinstance(node, dict):
        return
    for key in ("scope", "rss_measurement", "cpu_scope"):
        if isinstance(node.get(key), str):
            context[key] = node[key]
    if "full_inference" in node:
        if type(node["full_inference"]) is not bool:
            raise ValueError("Invalid full-inference flag")
        context["full_inference"] = node["full_inference"]
    if isinstance(node.get("status"), str):
        context["measurement_status"] = node["status"]
    if "returncode" in node:
        if type(node["returncode"]) is not int:
            raise ValueError("Invalid command return code")
        context["command_returncode"] = node["returncode"]
    budget = node.get("metadata", {}).get("cpu_budget") if isinstance(node.get("metadata"), dict) else None
    if budget is not None:
        context["requested_cpu_budget"] = budget
    if node.get("requested_cpus") is not None:
        context["requested_cpu_budget"] = node["requested_cpus"]
    scheduler = node.get("scheduler", {})
    if isinstance(scheduler, dict):
        if scheduler.get("AllocCPUS") is not None:
            value = scheduler["AllocCPUS"]
            if isinstance(value, bool) or not str(value).isdigit() or int(value) <= 0:
                raise ValueError("Invalid allocated CPU count")
            context["allocated_cpus"] = int(value)
        if scheduler.get("JobIDRaw") is not None:
            context["job_id"] = scheduler["JobIDRaw"]
        for source, target in (("State", "scheduler_state"), ("ExitCode", "scheduler_exit_code")):
            if scheduler.get(source) is not None:
                context[target] = scheduler[source]
    measured = node.get("measurement", node)
    if isinstance(measured, dict):
        values = dict(wall_seconds=first(measured, ("elapsed_seconds", "wall_s")),
            user_cpu_seconds=first(measured, ("user_seconds", "user_cpu_s")),
            system_cpu_seconds=first(measured, ("system_seconds", "system_cpu_s")),
            total_cpu_seconds_reported=first(measured, ("cpu_s",)))
        peak, unit, scope = memory({**measured, **({"memory": node["memory"]} if "memory" in node else {})}, context)
        if "memory_scope" in node and "memory" not in node:
            scope = node["memory_scope"]
        if any(value is not None for value in values.values()) or peak is not None:
            budget = context.get("requested_cpu_budget")
            if budget is not None and (type(budget) is not int or budget <= 0):
                raise ValueError("Invalid requested CPU budget")
            yield dict({k: context[k] for k in COLUMNS if k in context}, **values,
                peak_memory_value=peak, peak_memory_unit=unit, memory_scope=scope,
                requested_cpu_budget=budget, allocated_cpus=context.get("allocated_cpus"),
                gnu_exit_status=measured.get("exit_status"), source_location=list(location),
                measurement_scope=context.get("scope", "As recorded at source location; no full-run scope inferred"))
    for key, value in node.items():
        if key in SCALAR_WALL:
            if value is not None:
                yield dict({k: context[k] for k in COLUMNS if k in context},
                    wall_seconds=number(value), source_location=[*location, key],
                    measurement_scope="Named scalar interval; may overlap subintervals; do not sum")
        elif key not in EXCLUDE and key != "measurement":
            yield from intervals(value, (*location, key), context)


def producer_metrics(method, cell, checked):
    if not method.startswith("full_search:") or cell.get("status") != "independent_content_readback_verified":
        return None
    pin = cell.get("source_snapshot")
    if not isinstance(pin, dict) or not {"path", "bytes", "sha256"} <= pin.keys():
        return None
    check(pin)
    checked.append(pin)
    producer = json.loads(Path(pin["path"]).read_text())
    if "arms" not in producer:
        return None
    if producer.get("dataset") != cell["dataset"]:
        raise ValueError("Full-search producer dataset does not match the admitted cell")
    arms = [row for row in producer.get("arms", []) if row.get("label") == method.split(":", 1)[1]]
    if not arms:
        raise ValueError("Missing full-search producer arm")
    if len(arms) != 1:
        raise ValueError("Ambiguous full-search producer arm")
    native = arms[0].get("native", {})
    metrics = native.get("metrics")
    if not isinstance(metrics, dict):
        return None
    if (native.get("products") != cell.get("products")
            or any(native.get(key) != value for key, value in cell["coverage"].items())
            or native.get("status") != "fresh_native_content_verified"
            or metrics.get("status") != "complete" or metrics.get("harness", {}).get("exit_code") != 0):
        raise ValueError("Full-search native metrics are not bound to the admitted cell")
    index = producer["arms"].index(arms[0])
    check(pin)
    return pin, metrics, ["arms", index, "native", "metrics"]


def referenced_metrics(cell, checked):
    resources = cell.get("resources")
    if not isinstance(resources, dict):
        return
    for key in ("native_metrics", "native_gnu_time", "scoring_recovery_gnu_time"):
        pin = resources.get(key)
        if not isinstance(pin, dict) or not {"path", "bytes", "sha256"} <= pin.keys():
            continue
        check(pin)
        checked.append(pin)
        if key == "native_metrics":
            metrics = json.loads(Path(pin["path"]).read_text())
            root_hogs = cell.get("products", {}).get("orthohmm_root_hogs.tsv")
            counts, coverage = metrics.get("counts", {}), cell.get("coverage", {})
            if (metrics.get("status") != "complete" or root_hogs is None
                    or metrics.get("outputs", {}).get("root_hogs") != root_hogs
                    or counts.get("root_hogs") != coverage.get("groups")
                    or counts.get("ortholog_pairs") != coverage.get("pairs")):
                raise ValueError("Referenced native metrics are not bound to the admitted products/counts")
        else:
            # These explicit references use verbose GNU time; its existing reader rejects failed stages.
            metrics = verbose_time(Path(pin["path"]).read_text())
        check(pin)
        yield key, pin, metrics


def build(summary, source_path, checked):
    rows = []
    for method_index, method in enumerate(summary["methods"]):
        present = set()
        for index, cell in enumerate(method["cells"]):
            if cell["dataset"] not in DATASETS:
                raise ValueError("Unknown resource dataset")
            present.add(cell["dataset"])
            base = dict(method=method["key"], label=method["label"], dataset=cell["dataset"], cell_index=index,
                cell_status=cell["status"], cell_scope=cell.get("scope"), job_id=cell.get("producer_job_id"))
            location = ("methods", method_index, "cells", index, "resources")
            context = {"job_id": cell["producer_job_id"]} if cell.get("producer_job_id") is not None else {}
            original = list(intervals(cell.get("resources"), location, context))
            for value in original:
                rows.append({**base, "origin": "cell_resource", "source_document": str(source_path), **value})
            extra = producer_metrics(method["key"], cell, checked)
            added = []
            if extra is not None:
                pin, metrics, location = extra
                added = list(intervals(metrics, tuple(location), context))
                for value in added:
                    rows.append({**base, "origin": "bound_producer_native_metrics",
                                 "source_document": pin["path"], **value})
            for key, pin, metrics in referenced_metrics(cell, checked):
                values = list(intervals(metrics, context=context))
                added.extend(values)
                for value in values:
                    rows.append({**base, "origin": "explicit_resource_reference",
                        "source_document": pin["path"], "reference_location": ["methods", method_index,
                            "cells", index, "resources", key], **value,
                        "measurement_scope": cell.get("resources", {}).get("scope") or
                            "Explicit retained " + key + " reference; not inferred to be full-pipeline timing"})
            if not original and not added:
                rows.append({**base, "origin": "missing_resource", "source_document": str(source_path),
                    "source_location": ["methods", method_index, "cells", index],
                    "measurement_scope": "No retained resource interval for this historical cell; unknown is not zero"})
        for dataset in DATASETS:
            if dataset not in present:
                rows.append(dict(method=method["key"], label=method["label"], dataset=dataset,
                    cell_status="not_evaluated_no_cell", origin="missing_dataset_cell", source_document=str(source_path),
                    source_location=["methods", method_index],
                    measurement_scope="Dataset was not evaluated for this configuration; no resource value imputed"))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Preserve prior resource exports")
    source = record(args.comparison)
    summary = json.loads(args.comparison.read_text())
    checked = [source, record(__file__), *[record(Path(__file__).with_name(name)) for name in
        ("prepare_ob_candidate_neighborhood.py", "summarize_matched_resources.py", "gnu_time_companion.py")]]
    rows = build(summary, Path(source["path"]), checked)
    for pin in checked:
        check(pin)
    args.output.mkdir(parents=True)
    table = args.output / "resources.tsv"
    with table.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if row.get(key) is None else
                json.dumps(row[key]) if key in ("source_location", "reference_location") else row[key] for key in COLUMNS})
    result = dict(status="hillclimb_resource_intervals_exported_not_controlled_speed_ranking", source=source,
        rows=rows, table=record(table), checked=list({p["path"]: p for p in checked}.values()),
        limitations=["Intervals overlap and can be reused across configurations; rows are not independent runs and must not be summed.",
            "Native, cached, scoring, conversion, failed and recovery scopes remain separate; no isolated speed ranking on the shared host.",
            "Missing measurements or dataset cells remain NA, not zero. CPU totals are reported only when explicitly retained, not imputed from wall time.",
            "Memory units and process/parent/sampled-tree scopes are preserved. Exact byte values take precedence over duplicate rounded GiB fields.",
            "Source JSON locations and statuses identify each observation; inherited metadata is not new historical execution attestation."],
        defaults_promoted=False, goal_completed=False)
    with (args.output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"rows": len(rows), "methods": len(summary["methods"]), "table": str(table)}))


if __name__ == "__main__":
    main()

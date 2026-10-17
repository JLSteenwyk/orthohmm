"""Independent direct-location resource readback; no exporter implementation calls."""

import collections
import csv
import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET


def pin(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            digest.update(block)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def follow(value, location):
    for key in location:
        value = value[key]
    return value


def document(path):
    path = Path(path)
    if path.suffix == ".json":
        return json.loads(path.read_text())
    labels = {"User time (seconds)": "user_seconds", "System time (seconds)": "system_seconds",
        "Maximum resident set size (kbytes)": "max_process_rss_kib", "Exit status": "exit_status",
        "Elapsed (wall clock) time (h:mm:ss or m:ss)": "elapsed_seconds"}
    result = {}
    for line in path.read_text().splitlines():
        for label, key in labels.items():
            if line.strip().startswith(label + ": "):
                assert key not in result
                raw = line.strip()[len(label) + 2:]
                if key == "elapsed_seconds":
                    result[key] = sum(float(part) * 60 ** i for i, part in enumerate(reversed(raw.split(":"))))
                else:
                    result[key] = float(raw) if key in ("user_seconds", "system_seconds") else int(raw)
    assert len(result) == 5
    return result


def main(root):
    report_path = root / "report.json"
    report = json.loads(report_path.read_text())
    expected = {p["path"]: p for p in report["checked"]}
    expected[report["table"]["path"]] = report["table"]
    for path, item in expected.items():
        assert pin(path) == item, path
    data = {path: document(path) for path in expected if path.endswith((".json", ".txt"))}
    summary = data[report["source"]["path"]]
    methods = {method["key"]: method for method in summary["methods"]}
    row_counts = collections.Counter()
    memory_units = collections.Counter()
    origins = collections.Counter()
    missing = []
    scalar_fields = {"total_incremental_wall_s", "total_full_pipeline_and_scoring_wall_s", "total_pilot_wall_s",
        "native_and_scoring_wall_s", "candidate_wall_s", "scoring_recovery_wall_s"}
    aliases = {"wall_seconds": ("elapsed_seconds", "wall_s"),
        "user_cpu_seconds": ("user_seconds", "user_cpu_s"),
        "system_cpu_seconds": ("system_seconds", "system_cpu_s"), "total_cpu_seconds_reported": ("cpu_s",)}
    memory_keys = [("max_process_rss_kib", "KiB"), ("max_reported_process_rss_kib", "KiB"),
        ("peak_process_tree_rss_bytes", "bytes"), ("peak_process_tree_rss_gib", "GiB"),
        ("peak_process_rss_GiB", "GiB"), ("parent_peak_rss_GiB", "GiB")]
    for row in report["rows"]:
        method = methods[row["method"]]
        origins[row["origin"]] += 1
        if row["origin"] == "missing_dataset_cell":
            assert row["dataset"] not in {c["dataset"] for c in method["cells"]}
            assert row["cell_status"] == "not_evaluated_no_cell"
        else:
            cell = method["cells"][row["cell_index"]]
            row_counts[(row["method"], row["cell_index"])] += 1
            assert (cell["dataset"], cell["status"]) == (row["dataset"], row["cell_status"])
        if row["origin"].startswith("missing"):
            assert all(row.get(key) is None for key in [*aliases, "peak_memory_value"])
            missing.append({key: row.get(key) for key in ("method", "dataset", "cell_status", "origin")})
            continue
        source = data[row["source_document"]]
        location = row["source_location"]
        node = follow(source, location)
        scalar = bool(location and location[-1] in scalar_fields)
        measured = node.get("measurement", node) if isinstance(node, dict) else {}
        for target, names in aliases.items():
            observed = node if scalar and target == "wall_seconds" else None
            if not scalar:
                observed = next((measured[key] for key in names if key in measured), None)
            assert row.get(target) == observed, (row["method"], location, target)
        if not scalar:
            explicit = node.get("memory")
            if isinstance(explicit, dict):
                peak, unit, scope = explicit.get("value"), explicit.get("unit"), explicit.get("scope")
                assert row.get("memory_scope") == scope
            else:
                selected = next(((measured[key], unit) for key, unit in memory_keys if key in measured), (None, None))
                peak, unit = selected
                if "memory_scope" in node:
                    assert row.get("memory_scope") == node["memory_scope"]
            assert row.get("peak_memory_value") == peak
            assert row.get("peak_memory_unit") == unit
            assert row.get("gnu_exit_status") == measured.get("exit_status")
            memory_units[str(unit)] += 1
        start = {"cell_resource": 5, "bound_producer_native_metrics": 4,
                 "explicit_resource_reference": 0}[row["origin"]]
        ancestors = [follow(source, location[:i]) for i in range(start, len(location) + 1)]
        expected_job = cell.get("producer_job_id")
        budget = allocation = command_code = state = exit_code = full_inference = cpu_scope = None
        for ancestor in ancestors:
            if not isinstance(ancestor, dict):
                continue
            if isinstance(ancestor.get("metadata"), dict):
                budget = ancestor["metadata"].get("cpu_budget", budget)
            budget = ancestor.get("requested_cpus", budget)
            full_inference = ancestor.get("full_inference", full_inference)
            cpu_scope = ancestor.get("cpu_scope", cpu_scope)
            command_code = ancestor.get("returncode", command_code)
            scheduler = ancestor.get("scheduler", {})
            if not isinstance(scheduler, dict):
                continue
            expected_job = scheduler.get("JobIDRaw", expected_job)
            allocation = int(scheduler["AllocCPUS"]) if "AllocCPUS" in scheduler else allocation
            state = scheduler.get("State", state)
            exit_code = scheduler.get("ExitCode", exit_code)
        for key, value in (("job_id", expected_job), ("requested_cpu_budget", budget), ("allocated_cpus", allocation),
                ("command_returncode", command_code), ("scheduler_state", state), ("scheduler_exit_code", exit_code),
                ("full_inference", full_inference), ("cpu_scope", cpu_scope)):
            assert row.get(key) == value, (row["method"], location, key, row.get(key), value)
        if row["origin"] == "explicit_resource_reference":
            reference = follow(summary, row["reference_location"])
            assert reference == expected[row["source_document"]]
        elif row["origin"] == "bound_producer_native_metrics":
            assert cell["source_snapshot"] == expected[row["source_document"]]
            native = follow(source, location[:2])["native"]
            assert native["products"] == cell["products"]
            assert all(native[key] == value for key, value in cell["coverage"].items())
    cells = {(m["key"], i) for m in summary["methods"] for i, _ in enumerate(m["cells"])}
    assert cells == set(row_counts)
    assert len(cells) == 91 and len(methods) == 28
    with Path(report["table"]["path"]).open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        table = list(reader)
        columns = reader.fieldnames
    assert len(table) == len(report["rows"])
    for table_row, json_row in zip(table, report["rows"]):
        for key in columns:
            value = json_row.get(key)
            text = "NA" if value is None else json.dumps(value) if isinstance(value, list) else str(value)
            assert table_row[key] == text, key
    tests = Path("benchmarks/results/hillclimb_resource_export_tests_v2_20261016.xml")
    suites = list(ET.parse(tests).getroot().iter("testsuite"))
    counts = {key: sum(int(s.attrib.get(key, 0)) for s in suites) for key in ("tests", "failures", "errors", "skipped")}
    assert counts == dict(tests=63, failures=0, errors=0, skipped=0)
    for path, item in expected.items():
        assert pin(path) == item, path
    receipt = dict(status="independent_direct_source_and_tsv_resource_readback_verified", methods=len(methods),
        historical_cells=len(cells), rows=len(table), origins=dict(origins), memory_units=dict(memory_units),
        missing=missing, checked=list(expected.values()), report=pin(report_path), recipe=pin(__file__),
        tests={"counts": counts, "xml": pin(tests)}, defaults_promoted=False, goal_completed=False,
        limitations=["Direct source values, aliases, memory units, scopes, inherited CPU/job metadata and TSV serialization verified; no controlled timing comparison or new inference.",
            "Every historical cell is represented. Missing intervals and six unevaluated configuration/dataset combinations are explicit, not estimated.",
            "The two empty large-copy-component OrthoBench/ThreeKingdoms resource records correspond to inactive-policy native-output reuse, not new native inference.",
            "No costs are summed. Parent, process and sampled-tree memory remain distinct. Source identity does not prove historical input consumption."])
    with (root / "independent_readback.json").open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"rows": len(table), "historical_cells": len(cells), "origins": dict(origins), "tests": counts}))


if __name__ == "__main__":
    main(Path(sys.argv[1]))

import copy
import csv
import json
import sys

import pytest

from benchmark_tools import export_hillclimb_resources as export
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def summary(resources=None, **fields):
    cell = dict(dataset="OrthoBench", status="failed", resources=resources, **fields)
    return {"methods": [{"key": "test:arm", "label": "Test arm", "cells": [cell]}]}


@pytest.mark.parametrize("bad", [True, False, "10", -1, float("nan"), float("inf"), -float("inf")])
def test_invalid_measurements_rejected(bad):
    with pytest.raises(ValueError, match="finite nonnegative"):
        list(export.intervals({"wall_s": bad}))


def test_zero_is_a_measurement_but_none_is_unknown():
    rows = list(export.intervals({"wall_s": 0, "user_cpu_s": None}))
    assert len(rows) == 1
    assert rows[0]["wall_seconds"] == 0
    assert rows[0]["user_cpu_seconds"] is None
    assert list(export.intervals({"wall_s": None})) == []


def test_job_budget_scope_and_failure_state_inherited_without_summing():
    resources = {"scope": "cached replay", "metadata": {"cpu_budget": 16},
        "scheduler": {"AllocCPUS": "32", "State": "FAILED", "ExitCode": "1:0", "JobIDRaw": "123"},
        "native": {"returncode": -11, "resources": {"elapsed_seconds": 10, "exit_status": 139,
            "user_seconds": 12, "system_seconds": 2, "max_process_rss_kib": 100}},
        "total_incremental_wall_s": 11}
    rows = list(export.intervals(resources))
    assert len(rows) == 2
    native = next(r for r in rows if "gnu_exit_status" in r)
    assert native["source_location"] == ["native", "resources"]
    assert native["job_id"] == "123"
    assert native["requested_cpu_budget"] == 16
    assert native["allocated_cpus"] == 32
    assert native["scheduler_state"] == "FAILED"
    assert native["scheduler_exit_code"] == "1:0"
    assert native["command_returncode"] == -11
    assert native["gnu_exit_status"] == 139
    assert native["measurement_scope"] == "cached replay"
    assert native["total_cpu_seconds_reported"] is None
    assert next(r for r in rows if r["wall_seconds"] == 11)["job_id"] == "123"


def test_sibling_scope_and_job_do_not_leak():
    rows = list(export.intervals([
        {"scope": "conversion", "scheduler": {"JobIDRaw": "123"}, "wall_s": 2}, {"wall_s": 3}]))
    assert rows[0]["job_id"] == "123"
    assert "job_id" not in rows[1]
    assert rows[1]["measurement_scope"] != "conversion"


def test_cell_job_fallback_is_preserved_and_scheduler_job_takes_precedence():
    source = summary([{"measurement": {"elapsed_seconds": 2}},
        {"scheduler": {"JobIDRaw": "child"}, "measurement": {"elapsed_seconds": 3}}], producer_job_id="parent")
    rows = export.build(source, "summary.json", [])
    assert [r["job_id"] for r in rows[:2]] == ["parent", "child"]
    assert len(rows) == 4


@pytest.mark.parametrize("key,unit,scope", [
    ("max_process_rss_kib", "KiB", "process"),
    ("peak_process_tree_rss_bytes", "bytes", "tree"),
    ("peak_process_tree_rss_gib", "GiB", "tree"),
    ("peak_process_rss_GiB", "GiB", "process"),
    ("parent_peak_rss_GiB", "GiB", "parent"),
])
def test_memory_units_and_scopes_not_converted(key, unit, scope):
    row, = export.intervals({key: 100})
    assert row["peak_memory_value"] == 100
    assert row["peak_memory_unit"] == unit
    assert scope in row["memory_scope"].lower()


def test_exact_bytes_precede_duplicate_rounded_gib_and_keep_sampling_scope():
    row, = export.intervals({"peak_process_tree_rss_bytes": 1073741824,
        "peak_process_tree_rss_gib": 1, "rss_measurement": "sampled_sum_of_linux_proc_tree_rss"})
    assert row["peak_memory_value"] == 1073741824
    assert row["peak_memory_unit"] == "bytes"
    assert row["memory_scope"] == "sampled_sum_of_linux_proc_tree_rss"


def test_explicit_parent_memory_and_cpu_scope_are_retained():
    row, = export.intervals({"full_inference": False, "scope": "driver only",
        "cpu_scope": "excludes Docker tasks", "measurement": {"elapsed_seconds": 3},
        "memory": {"value": 99, "unit": "KiB", "scope": "parent only"}})
    assert row["source_location"] == []
    assert row["peak_memory_value"] == 99
    assert row["memory_scope"] == "parent only"
    assert row["full_inference"] is False
    assert row["cpu_scope"] == "excludes Docker tasks"


def test_outer_memory_scope_is_not_lost_when_measurement_is_nested():
    row, = export.intervals({"memory_scope": "not aggregate", "measurement": {"max_process_rss_kib": 10}})
    assert row["memory_scope"] == "not aggregate"


@pytest.mark.parametrize("field,value", [
    ("requested_cpus", True), ("requested_cpus", 0), ("requested_cpus", "16"),
    ("full_inference", "yes"), ("returncode", "1"),
])
def test_invalid_execution_metadata_rejected(field, value):
    with pytest.raises(ValueError):
        list(export.intervals({field: value, "wall_s": 2}))


@pytest.mark.parametrize("value", [True, 0, "0", -1, "1.0", "not known"])
def test_invalid_scheduler_allocation_rejected(value):
    with pytest.raises(ValueError, match="allocated CPU"):
        list(export.intervals({"scheduler": {"AllocCPUS": value}, "wall_s": 2}))


def test_metadata_and_retained_measurement_duplicates_are_not_measured_twice():
    rows = list(export.intervals({"wall_s": 1, "retained_measurements": {"wall_s": 1},
        "metadata": {"wall_s": 50}, "host_before": {"wall_s": 999},
        "stages": {"alignment": {"wall_s": 0.5}}}))
    assert [r["wall_seconds"] for r in rows] == [1, 0.5]


def test_failed_and_not_evaluated_cells_stay_explicit_not_zero():
    rows = export.build(summary(), "summary.json", [])
    assert len(rows) == 3
    assert rows[0]["cell_status"] == "failed"
    assert rows[0]["origin"] == "missing_resource"
    assert rows[1]["cell_status"] == rows[2]["cell_status"] == "not_evaluated_no_cell"
    assert all("wall_seconds" not in row for row in rows)


def test_unknown_dataset_rejected():
    source = summary()
    source["methods"][0]["cells"][0]["dataset"] = "other"
    with pytest.raises(ValueError, match="Unknown resource dataset"):
        export.build(source, "summary.json", [])


def producer(tmp_path):
    products = {"native": {"path": "product", "bytes": 5, "sha256": "identity"}}
    coverage = {"genes": 4, "groups": 2, "pairs": 3, "species": 2}
    metrics = {"wall_s": 10, "metadata": {"cpu_budget": 16},
        "status": "complete", "harness": {"exit_code": 0}}
    native = dict(products=products, **coverage, status="fresh_native_content_verified", metrics=metrics)
    path = tmp_path / "producer.json"
    path.write_text(json.dumps({"dataset": "OrthoBench", "arms": [{"label": "arm", "native": native}]}))
    cell = dict(dataset="OrthoBench", status="independent_content_readback_verified",
        source_snapshot=record(path), coverage=coverage, products=products,
        resources={"total_full_pipeline_and_scoring_wall_s": 12})
    return cell, path


def test_bound_producer_metrics_supplement_but_do_not_replace_scalar(tmp_path):
    cell, path = producer(tmp_path)
    source = {"methods": [{"key": "full_search:arm", "label": "Arm", "cells": [cell]}]}
    checked = []
    rows = export.build(source, "summary.json", checked)
    assert len(rows) == 4
    assert rows[0]["wall_seconds"] == 12
    assert rows[1]["wall_seconds"] == 10
    assert rows[1]["source_document"] == str(path)
    assert rows[1]["source_location"] == ["arms", 0, "native", "metrics"]
    assert rows[1]["requested_cpu_budget"] == 16
    assert checked == [record(path)]


@pytest.mark.parametrize("field", ["products", "coverage", "dataset"])
def test_mismatched_producer_cannot_supply_another_arms_measurements(tmp_path, field):
    cell, _ = producer(tmp_path)
    cell[field] = "wrong" if field == "dataset" else {"wrong": 1}
    with pytest.raises(ValueError):
        export.producer_metrics("full_search:arm", cell, [])


def test_changed_pinned_producer_rejected(tmp_path):
    cell, path = producer(tmp_path)
    path.write_text("{}")
    with pytest.raises(ValueError, match="identity changed"):
        export.producer_metrics("full_search:arm", cell, [])


def test_ambiguous_producer_arm_rejected(tmp_path):
    cell, path = producer(tmp_path)
    doc = json.loads(path.read_text())
    doc["arms"].append(copy.deepcopy(doc["arms"][0]))
    path.write_text(json.dumps(doc))
    cell["source_snapshot"] = record(path)
    with pytest.raises(ValueError, match="Ambiguous"):
        export.producer_metrics("full_search:arm", cell, [])


def test_only_full_search_verified_cells_use_producer_supplement(tmp_path):
    cell, _ = producer(tmp_path)
    assert export.producer_metrics("other:arm", cell, []) is None
    cell["status"] = "failed"
    assert export.producer_metrics("full_search:arm", cell, []) is None


def test_recovery_receipt_is_not_mistaken_for_an_arm_producer(tmp_path):
    cell, path = producer(tmp_path)
    path.write_text(json.dumps({"native_summary": {}, "coverage": {}}))
    cell["source_snapshot"] = record(path)
    assert export.producer_metrics("full_search:arm", cell, []) is None


def test_missing_arm_in_an_arm_producer_is_rejected(tmp_path):
    cell, _ = producer(tmp_path)
    with pytest.raises(ValueError, match="Missing full-search producer arm"):
        export.producer_metrics("full_search:another", cell, [])


def test_explicit_metrics_reference_requires_matching_products_and_counts(tmp_path):
    root_hogs = {"path": "hogs", "bytes": 10, "sha256": "identity"}
    path = tmp_path / "metrics.json"
    path.write_text(json.dumps({"status": "complete", "outputs": {"root_hogs": root_hogs},
        "counts": {"root_hogs": 2, "ortholog_pairs": 3}, "wall_s": 10}))
    cell = summary({"native_metrics": record(path)}, coverage={"groups": 2, "pairs": 3},
        products={"orthohmm_root_hogs.tsv": root_hogs})["methods"][0]["cells"][0]
    checked = []
    key, pin, metrics = next(export.referenced_metrics(cell, checked))
    assert (key, pin, metrics["wall_s"]) == ("native_metrics", record(path), 10)
    assert checked == [record(path)]
    cell["coverage"]["pairs"] = 4
    with pytest.raises(ValueError, match="not bound"):
        list(export.referenced_metrics(cell, []))


def test_cli_preserves_unknown_as_na_and_refuses_to_overwrite(tmp_path, monkeypatch):
    source = tmp_path / "summary.json"
    source.write_text(json.dumps(summary({"wall_s": 0}, producer_job_id="123")))
    output = tmp_path / "export"
    monkeypatch.setattr(sys, "argv", ["export", "--comparison", str(source), "--output", str(output)])
    export.main()
    report = json.loads((output / "report.json").read_text())
    assert report["source"] == record(source)
    assert report["table"] == record(output / "resources.tsv")
    assert not report["goal_completed"] and not report["defaults_promoted"]
    with (output / "resources.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 3
    assert rows[0]["wall_seconds"] == "0"
    assert rows[0]["user_cpu_seconds"] == rows[0]["peak_memory_value"] == "NA"
    assert rows[0]["job_id"] == "123"
    with pytest.raises(FileExistsError):
        export.main()

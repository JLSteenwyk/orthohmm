"""Keep cached-stage evidence separate from absent full-pipeline measurements."""

import copy
import json
import hashlib
from pathlib import Path

import pytest

from benchmark_tools import collect_factorial_resources as costs


def fixture():
    plan = {
        "cells": [{"label": label, "profile_expansion": label[1] == "1",
                   "candidate_expansion": label[4] == "1", "reconciliation": label[7] == "1",
                   "runtime_kind": "incremental_cached_replay"} for label in costs.CELLS],
        "candidate_arms": {arm: {"incremental_preparation_seconds": 2.0} for arm in costs.ARMS},
    }
    metrics = {"status": "complete", "metadata": {"cpu_budget": 32},
               "rss_measurement": "sampled_sum_of_linux_proc_tree_rss",
               "wall_s": 10.0, "user_cpu_s": 100.0, "system_cpu_s": 20.0,
               "peak_process_tree_rss_bytes": 2**30}
    observed = {label: copy.deepcopy(metrics) for label in costs.CELLS if label.endswith("r1")}
    batches = {label: {"State": "FAILED", "ExitCode": "1:0"} for label in observed}
    return plan, observed, batches


def test_no_imputed_full_costs_shared_preparation_not_doubled_and_failures_remain():
    plan, measured, batches = fixture()
    rows = costs.project("OrthoBench", plan, measured, batches)
    assert [row["cell"] for row in rows] == list(costs.CELLS)
    assert all(row["full_pipeline_wall_s"] is row["full_pipeline_cpu_s"]
               is row["full_pipeline_peak_memory_bytes"] is None for row in rows)
    assert all(row["candidate_preparation_cpu_s"] is row["candidate_preparation_peak_bytes"] is None
               for row in rows)
    for index in range(0, 8, 2):
        off, on = rows[index:index + 2]
        assert off["candidate_preparation_arm"] == on["candidate_preparation_arm"]
        assert off["candidate_preparation_wall_s"] == on["candidate_preparation_wall_s"] == 2.0
        assert off["reconciliation_measurement_status"] == "not_applicable"
        assert off["reconciliation_wall_s"] is off["original_batch"] is None
        assert on["original_batch"] == {"State": "FAILED", "ExitCode": "1:0"}
        assert on["reconciliation_mean_cpu_cores"] == 12.0
    rendered = costs.render({"rows": rows, "limitations": ["Do not infer isolated speed."]})
    assert rendered.count("| NA |") >= 8 and "Do not infer isolated speed." in rendered


@pytest.mark.parametrize("damage", ["missing_cell", "duplicate_cell", "flags", "flag_type",
                                   "scope", "arm", "preparation", "missing_metrics", "extra_metrics", "batch"])
def test_bad_factorial_inventory_or_scope_is_rejected(damage):
    plan, measured, batches = fixture()
    if damage == "missing_cell":
        plan["cells"].pop()
    elif damage == "duplicate_cell":
        plan["cells"][-1] = plan["cells"][0]
    elif damage == "flags":
        plan["cells"][0]["reconciliation"] = True
    elif damage == "flag_type":
        plan["cells"][0]["profile_expansion"] = 0
    elif damage == "scope":
        plan["cells"][0]["runtime_kind"] = "full_pipeline"
    elif damage == "arm":
        plan["candidate_arms"].pop("p0_c0")
    elif damage == "preparation":
        plan["candidate_arms"]["p0_c0"]["incremental_preparation_seconds"] = 0
    elif damage == "missing_metrics":
        measured.pop("p0_c0_r1")
    elif damage == "extra_metrics":
        measured["p0_c0_r0"] = measured["p0_c0_r1"]
    else:
        batches.pop("p0_c0_r1")
    with pytest.raises(ValueError):
        costs.project("OrthoBench", plan, measured, batches)


@pytest.mark.parametrize("key,value", [
    ("wall_s", 0), ("wall_s", True), ("wall_s", float("inf")), ("wall_s", float("nan")),
    ("user_cpu_s", -1), ("system_cpu_s", "20"), ("peak_process_tree_rss_bytes", 1.5),
    ("peak_process_tree_rss_bytes", False), ("status", "failed"),
    ("rss_measurement", "cgroup_lifetime_peak"),
])
def test_invalid_resource_observations_are_not_relabelled(key, value):
    _, measured, _ = fixture()
    metrics = measured["p0_c0_r1"]
    metrics[key] = value
    with pytest.raises(ValueError):
        costs.metrics_values(metrics)


def test_record_binding_and_repository_boundary(tmp_path):
    path = tmp_path / "metrics.json"
    path.write_text(json.dumps({"wall_s": 10.0}))
    ref = costs.record(tmp_path, path)
    assert costs.load(tmp_path, ref)[0] == {"wall_s": 10.0}
    path.write_text(json.dumps({"wall_s": 11.0}))
    with pytest.raises(ValueError, match="checksum"):
        costs.load(tmp_path, ref)
    with pytest.raises(ValueError, match="outside"):
        costs.record(tmp_path, tmp_path.parent)


def test_existing_output_refused_without_reading_or_writing_inputs(tmp_path):
    marker = tmp_path / "marker"
    marker.write_text("preserve")
    with pytest.raises(FileExistsError):
        costs.write(tmp_path / "nonexistent_inputs", tmp_path)
    assert marker.read_text() == "preserve"
    assert list(tmp_path.iterdir()) == [marker]


def test_actual_sixteen_cell_export_agrees_with_independent_direct_readback():
    root = Path(__file__).resolve().parents[2]
    base = root / "benchmark_tools/results/factorial_retained_resources_20261004"
    report = json.loads((base / "resources.json").read_text())
    assert report["cells"] == 16 and report["recorded_reconciliation_cells"] == 8
    assert report["full_pipeline_cost_cells_available"] == 0
    assert len(report["inputs"]) == 15
    for ref in [report["source"], *report["inputs"]]:
        data = (root / ref["path"]).read_bytes()
        assert (len(data), hashlib.sha256(data).hexdigest()) == (ref["bytes"], ref["sha256"])
    ob = json.loads((root / costs.FIXED_INPUTS["ob_results"][0]).read_text())
    plans = {"OrthoBench": json.loads((root / costs.FIXED_INPUTS["ob_preparation"][0]).read_text()),
             "Corrected QfO": json.loads((root / costs.FIXED_INPUTS["qfo_preparation"][0]).read_text())}
    rows = {(r["dataset"], r["cell"]): r for r in report["rows"]}
    assert set(rows) == {(dataset, cell) for dataset in plans for cell in costs.CELLS}
    for (dataset, cell), row in rows.items():
        arm = plans[dataset]["candidate_arms"][cell[:-3]]
        assert row["candidate_preparation_wall_s"] == arm["incremental_preparation_seconds"]
        assert row["full_pipeline_wall_s"] is row["full_pipeline_cpu_s"] is row["full_pipeline_peak_memory_bytes"] is None
        if cell.endswith("r0"):
            assert row["reconciliation_wall_s"] is row["reconciliation_peak_sampled_tree_rss_bytes"] is None
            continue
        if dataset == "OrthoBench":
            native = ob["native_validation"][cell]
        else:
            job = costs.QFO_ADMISSIONS[cell][0]
            admission = json.loads((root / ("benchmark_tools/results/qfo_corrected_factorial_native_admission_%d.json" % job)).read_text())
            native = admission["native_group_integrity"]
        metrics = json.loads(Path(native["native_metrics"]["path"]).read_text())
        assert row["reconciliation_wall_s"] == metrics["wall_s"]
        assert row["reconciliation_user_cpu_s"] == metrics["user_cpu_s"]
        assert row["reconciliation_system_cpu_s"] == metrics["system_cpu_s"]
        assert row["reconciliation_peak_sampled_tree_rss_bytes"] == metrics["peak_process_tree_rss_bytes"]
        assert row["reconciliation_mean_cpu_cores"] == (metrics["user_cpu_s"] + metrics["system_cpu_s"]) / metrics["wall_s"]
        assert row["original_batch"] == native["integrity"]["scheduler"]

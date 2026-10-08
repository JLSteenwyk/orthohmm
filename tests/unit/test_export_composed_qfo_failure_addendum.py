"""Failure reporting arithmetic and no-admission semantics, not live scores."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import export_composed_qfo_failure_addendum as current
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture():
    failed = dict(name="fas_benchmark (1)", status="FAILED", exit="137")
    failure = dict(schema="native11_composed_assessment_failure_readback_v1",
        status="failed_scoring_retained_not_admitted", job_id=24038, native_job_id=23985,
        native_index=11, cell="p1_c1_r0", output_inventory_reproduced=True,
        accuracy_admitted=False, admission_submitted=False, automatic_retry=False,
        native_inference_reexecuted=False, publication_ready=False, execution_status="failed", execution_exit_code=1,
        failed_tasks=[failed], trace_tasks=[failed, *[dict(name=name, status="COMPLETED", exit="0") for name in (
            "vgnc_benchmark (1)", "reference_genetrees_benchmark (SwissTrees)",
            "reference_genetrees_benchmark (TreeFam-A)", "ec_benchmark (1)", "go_benchmark (1)")]],
        accounting=dict(rows=[dict(JobIDRaw="24038", State="OUT_OF_MEMORY", ExitCode="0:125", ReqMem="32G",
            AllocCPUS="8", Elapsed="00:16:06", MaxRSS="", MaxVMSize="")]))
    review = dict(schema="native11_composed_terminal_review_v1", status="native_success", index=11,
        job_id=23985, cell="p1_c1_r0", terminal_reviewed=True, composed_full_review_complete=True,
        native_outputs_validated=True, primary_resources_replayed=True, shared_host_resources_reviewed=True,
        resource_scopes=current.SCOPES, accuracy_evaluated=False, uncontended_timing=False,
        resources=dict(wall_seconds=60568.563058126, cpu_seconds=1607563.508002, peak_memory_bytes=20407123968))
    conversion = dict(schema="native11_composed_qfo_conversion_v1", status="native11_composed_qfo_pairs_prepared_unscored",
        native_index=11, native_job_id=23985, job_id="24035", cell="p1_c1_r0", conversion_kind="group",
        semantics="cross-species group-derived clique pairs", accuracy_evaluated=False,
        native_inference_reexecuted=False, removed_mapping_pairs=0,
        retained_pairs=11755521, total_pairs=11755521, expected_pairs=11755521,
        pair_coverage=dict(input_accessions=984137, accessions_in_any_pair=585180, pair_rows=11755521,
            fraction_inputs_in_any_pair=585180 / 984137))
    return dict(failure=failure, review=review, conversion=conversion)


def test_failed_scoring_has_no_admitted_metric_or_mean():
    row = current.validate_reports(**fixture())
    assert row["accuracy_admitted"] is False
    assert row["secondary_six_metric_mean"] is None
    assert row["admitted_endpoint_score_count"] == 0
    assert row["completed_endpoint_task_count"] == 5
    assert row["scoring_peak_memory_bytes"] is None
    assert row["relation_coverage"] == 585180 / 984137
    assert row["resources"]["peak_memory_bytes"] != row["scoring_memory_limit_bytes"]


@pytest.mark.parametrize("component,field,value", [
    ("failure", "schema", "admitted"), ("failure", "accuracy_admitted", True),
    ("failure", "admission_submitted", True), ("failure", "automatic_retry", True),
    ("failure", "output_inventory_reproduced", False), ("failure", "failed_tasks", []),
    ("failure", "execution_status", "success"), ("failure", "native_index", 12),
    ("review", "status", "failure"), ("review", "native_outputs_validated", False),
    ("review", "primary_resources_replayed", False), ("review", "resource_scopes", {}),
    ("review", "uncontended_timing", True), ("review", "accuracy_evaluated", True),
    ("review", "resources", dict(wall_seconds=1, cpu_seconds=1, peak_memory_bytes=1.0)),
    ("conversion", "conversion_kind", "native"), ("conversion", "accuracy_evaluated", True),
    ("conversion", "removed_mapping_pairs", 1), ("conversion", "expected_pairs", 1)])
def test_inconsistent_failure_coverage_or_scope_refuses(component, field, value):
    values = fixture()
    values[component][field] = value
    with pytest.raises(ValueError):
        current.validate_reports(**values)


def test_zero_imputation_or_partial_trace_cannot_claim_completion():
    values = fixture()
    values["failure"]["trace_tasks"] = values["failure"]["trace_tasks"][:-1]
    with pytest.raises(ValueError):
        current.validate_reports(**values)


def test_machine_readable_and_markdown_export_are_separate_and_unscored(tmp_path, monkeypatch):
    values = fixture()
    manuscript = tmp_path / "frozen_manuscript.md"
    manuscript.write_text("frozen manuscript fixture\n")
    refs = {"manuscript": record(manuscript)}
    for name, value in values.items():
        path = tmp_path / (name + ".json")
        path.write_text(json.dumps(value))
        refs[name] = record(path)
    monkeypatch.setattr(current, "ROOT", tmp_path)
    monkeypatch.setattr(current, "inputs", lambda: (deepcopy(refs), deepcopy(values)))
    output = tmp_path / "addendum"
    result = current.export(output)
    report = json.loads(Path(result["path"]).read_text())
    assert report["schema"] == "composed_native_qfo_failure_reporting_v1"
    assert report["new_scientific_admission"] is False
    assert report["frozen_manuscript_replaced"] is False
    assert report["publication_ready"] is False
    assert report["row"]["secondary_six_metric_mean"] is None
    assert "Unavailable\tUnavailable" in (output / "status.tsv").read_text()
    assert "coverage, not accuracy" in (output / "addendum.md").read_text()
    assert "exact peak RSS" in (output / "addendum.md").read_text()
    assert record(manuscript) == refs["manuscript"]
    with pytest.raises(ValueError, match="new direct"):
        current.export(output)

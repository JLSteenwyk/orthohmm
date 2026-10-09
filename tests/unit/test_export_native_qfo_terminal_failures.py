"""Terminal failures remain missing; reporting does not change admitted scores."""

from copy import deepcopy
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import export_native_qfo_terminal_failures as current
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture():
    rows = []
    for index, cell in zip(range(6, 13), current.CELLS):
        admitted = index in current.ADMITTED
        rows.append(dict(index=index, cell=cell, accuracy_admitted=admitted, status="frozen_snapshot_status",
            scores={name: 0.5 if admitted else None for name in current.ENDPOINTS},
            secondary_mean=0.5 if admitted else None, resources=None,
            timing_eligible=False if index == 7 else None))
    snapshot = dict(schema="allocated_native_qfo_scientific_reporting_snapshot_v1",
        supplied_admissions=4, publication_ready=False, rows=rows)
    resources = dict(wall_seconds=64526.0, cpu_seconds=1700000.0, peak_memory_bytes=22000000000)
    failed = dict(native_index=11, cell="p1_c1_r0", native_job_id=23985, conversion_job_id=24035,
        assessment_job_id=24038, inference_status="successful_composed_terminal_review",
        conversion_status="successful_unscored_group_cliques", scoring_status="OUT_OF_MEMORY",
        failed_endpoint="FAS", failed_endpoint_exit_code=137, completed_endpoint_task_count=5,
        admitted_endpoint_score_count=0, accuracy_admitted=False, secondary_six_metric_mean=None,
        submitted_pairs=11755521, input_accessions=984137, relation_accessions=585180,
        relation_coverage=585180 / 984137, prediction_semantics="cross-species group-derived clique pairs",
        resources=resources, resource_scopes=current.SCOPES)
    native11 = dict(schema="composed_native_qfo_failure_reporting_v1", row=failed,
        publication_ready=False, new_scientific_admission=False, automatic_retry=False,
        native_inference_reexecuted=False, frozen_manuscript_replaced=False)
    review = dict(schema="native12_composed_terminal_review_v1", status="native_failure_retained", index=12,
        job_id=24036, cell="p1_c1_r1", scheduler_state="FAILED", scheduler_exit_code="1:0",
        resource_scopes=current.SCOPES, resources=resources, terminal_reviewed=True,
        primary_resources_replayed=True, shared_host_resources_reviewed=True, native_outputs_validated=False,
        accuracy_evaluated=False, automatic_retry=False, publication_ready=False,
        uncontended_timing=False, scientific_timings_admitted=False)
    components = dict(scheduler=dict(verified=dict(fields=dict(JobState="FAILED", ExitCode="1:0"))),
        outputs_or_failure=dict(status="native_failure_retained", native_outcome="exited_nonzero",
            native_exit_code=-11, native_outputs_validated=False, accuracy_evaluated=False, automatic_retry=False),
        resource_replay=dict(schema="native12_composed_resource_replay_summary_v1",
            native_outcome="exited_nonzero", native_exit_code=-11, full_replay_executed=True,
            measured_matches_retained_wrapper=True),
        resources=dict(primary=resources, primary_scopes=current.SCOPES, native_outcome="exited_nonzero",
            native_exit_code=-11, shared_host_observation=True, uncontended_timing=False),
        environment=dict(status="shared_environment_replayed", sampled_environment_evidence_valid=True,
            uncontended_timing=False, background_cpu_used_for_eligibility=False,
            pressure_thresholds_used_for_eligibility=False),
        native_done=dict(exit_code=-11, timed_out=False),
        native_receipt=dict(schema="allocated_native_factorial_execution_v1", status="native_factorial_running",
            index=12, cell="p1_c1_r1"),
        missing_final_outputs={"metrics.json": True, "native/orthohmm_orthogroups.txt": True})
    return dict(snapshot=snapshot, native11=native11, review=review, components=components)


def test_preserves_all_four_admitted_rows_and_all_missing_scores():
    values = fixture()
    original = deepcopy(values)
    rows = current.validate(**values)
    assert values == original
    assert [row for row in rows if row["accuracy_admitted"]] == [
        row for row in original["snapshot"]["rows"] if row["accuracy_admitted"]]
    assert rows[3] == original["snapshot"]["rows"][3]
    assert sum(row["accuracy_admitted"] for row in rows) == 4
    assert all(value is None for row in rows if not row["accuracy_admitted"] for value in row["scores"].values())
    assert rows[5]["relation_coverage"] == 585180 / 984137
    assert rows[6]["resources"] == values["review"]["resources"]
    assert rows[6]["submitted_pairs"] is rows[6]["relation_coverage"] is None
    assert rows[6]["scoring_status"] == "not_run_failed_inference"


@pytest.mark.parametrize("component,field,value", [
    ("snapshot", "supplied_admissions", 5), ("snapshot", "publication_ready", True),
    ("native11", "new_scientific_admission", True), ("native11", "automatic_retry", True),
    ("review", "status", "native_success"), ("review", "terminal_reviewed", False),
    ("review", "primary_resources_replayed", False), ("review", "shared_host_resources_reviewed", False),
    ("review", "native_outputs_validated", True), ("review", "accuracy_evaluated", True),
    ("review", "scheduler_state", "COMPLETED"), ("review", "scheduler_exit_code", "0:0"),
    ("review", "job_id", 24080), ("review", "index", 11), ("review", "cell", "p1_c1_r0"),
    ("review", "resource_scopes", {}), ("review", "uncontended_timing", True),
    ("review", "scientific_timings_admitted", True),
])
def test_rejects_incompatible_failure_reporting(component, field, value):
    values = fixture()
    values[component][field] = value
    with pytest.raises(ValueError):
        current.validate(**values)


@pytest.mark.parametrize("component,field,value", [
    ("outputs_or_failure", "native_exit_code", 0), ("outputs_or_failure", "accuracy_evaluated", True),
    ("outputs_or_failure", "native_outcome", "timed_out"),
    ("resource_replay", "native_exit_code", 139), ("resource_replay", "full_replay_executed", False),
    ("resource_replay", "measured_matches_retained_wrapper", False),
    ("resources", "native_exit_code", 0), ("resources", "primary_scopes", {}),
    ("environment", "sampled_environment_evidence_valid", False),
    ("environment", "background_cpu_used_for_eligibility", True),
    ("environment", "pressure_thresholds_used_for_eligibility", True),
    ("native_done", "exit_code", 0), ("native_done", "timed_out", True),
    ("native_receipt", "status", "native_factorial_completed_pending_output_review"),
    ("missing_final_outputs", "metrics.json", False),
])
def test_rejects_wrapper_success_and_changed_contention_policy(component, field, value):
    values = fixture()
    values["components"][component][field] = value
    with pytest.raises(ValueError):
        current.validate(**values)


@pytest.mark.parametrize("field,value", [
    ("scoring_status", "COMPLETED"), ("admitted_endpoint_score_count", 5),
    ("completed_endpoint_task_count", 6), ("secondary_six_metric_mean", 0),
    ("relation_coverage", 0), ("resources", dict(wall_seconds=1, cpu_seconds=1, peak_memory_bytes=1.0)),
])
def test_rejects_native11_partial_score_admission_or_bad_coverage(field, value):
    values = fixture()
    values["native11"]["row"][field] = value
    with pytest.raises(ValueError):
        current.validate(**values)


@pytest.mark.parametrize("value", [0, float("nan"), float("inf"), True])
def test_resource_values_are_finite_positive_measurements(value):
    values = fixture()
    values["review"]["resources"]["wall_seconds"] = value
    with pytest.raises(ValueError):
        current.validate(**values)


@pytest.mark.parametrize("index", [9, 11, 12])
def test_missing_score_and_mean_imputation_refused(index):
    values = fixture()
    values["snapshot"]["rows"][index - 6]["scores"]["FAS"] = 0.0
    with pytest.raises(ValueError, match="imputed"):
        current.validate(**values)


@pytest.mark.parametrize("stdout", [
    "24080|RUNNING|0:0|2|128G\n", "24080|FAILED|1:0|2|128G\n",
    "24080|COMPLETED|0:0|2|32G\n", "24036|COMPLETED|0:0|2|128G\n", "",
])
def test_independent_reviewer_must_be_successful_terminal_job(stdout, monkeypatch):
    monkeypatch.setattr(current.subprocess, "run", lambda *args, **kwargs:
        SimpleNamespace(stdout=stdout, stderr="", returncode=0))
    with pytest.raises(ValueError, match="reviewer"):
        current.reviewer_accounting()


def test_reviewer_accounting_records_actual_observation(monkeypatch):
    observed = []

    def run(command, **kwargs):
        observed.append((command, kwargs))
        return SimpleNamespace(stdout="24080|COMPLETED|0:0|2|128G\n", stderr="", returncode=0)

    monkeypatch.setattr(current.subprocess, "run", run)
    result = current.reviewer_accounting()
    assert result["command"] == observed[0][0]
    assert observed[0][1]["check"] is True
    assert observed[0][1]["timeout"] == 5


def test_export_leaves_parent_unchanged_and_preserves_relative_link_location(tmp_path, monkeypatch):
    values = fixture()
    parent = tmp_path / "frozen.md"
    parent.write_text("old header\n## Abstract\nFour admitted cells.\n[Evidence](evidence.json)\n"
        + current.ANCHOR + "remainder unchanged\n")
    docs = dict(parent=parent.read_text(), snapshot=values["snapshot"], native11=values["native11"], review=values["review"])
    refs = {"parent": record(parent)}
    for key in ("snapshot", "native11", "review"):
        path = tmp_path / (key + ".json")
        path.write_text(json.dumps(docs[key]))
        refs[key] = record(path)
    monkeypatch.setattr(current, "ROOT", tmp_path)
    monkeypatch.setattr(current, "inputs", lambda *args: (deepcopy(docs), deepcopy(refs), deepcopy(values["components"])))
    monkeypatch.setattr(current, "reviewer_accounting", lambda: dict(fixture_only=True))
    destination = tmp_path / "reporting"
    result = current.export(Path(refs["review"]["path"]), refs["review"]["sha256"], destination)
    report = json.loads(Path(result["path"]).read_text())
    assert report["schema"] == "native_qfo_terminal_failure_reporting_v1"
    assert (report["admitted_endpoints"], report["missing_endpoints"]) == (24, 18)
    assert all(report[key] is False for key in ("publication_ready", "new_scientific_admission",
        "native_inference_reexecuted", "automatic_retry", "frozen_manuscript_replaced"))
    assert record(parent) == refs["parent"]
    manuscript = tmp_path / "reporting_manuscript.md"
    assert manuscript.exists() and manuscript.parent == parent.parent
    assert "[Evidence](evidence.json)" in manuscript.read_text()
    assert "remainder unchanged" in manuscript.read_text()
    with (destination / "scores.tsv").open() as stream:
        scores = list(csv.DictReader(stream, delimiter="\t"))
    assert len(scores) == 42
    assert sum(row["Value"] == "Unavailable" for row in scores) == 18
    assert "root cause and\ncrash location are unknown" in (destination / "addendum.md").read_text()
    with pytest.raises(ValueError, match="fresh direct"):
        current.export(Path(refs["review"]["path"]), refs["review"]["sha256"], destination)
    with pytest.raises(ValueError, match="fresh manuscript"):
        current.export(Path(refs["review"]["path"]), refs["review"]["sha256"], tmp_path / "other" / "reporting")

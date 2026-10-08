"""Explicit synthetic admissions with real score arithmetic, never live scores."""

from copy import deepcopy
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import export_composed_native_qfo_scientific_scores as current
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_native12_composed_review_binding import fixture as review_fixture
from tests.unit.test_native12_composed_qfo_assessment import (
    joined as assessment_joined, original_joined, short_root, execute, admit_fixture, put, runner)


@pytest.fixture
def joined(assessment_joined, monkeypatch):
    root, old_ref, stage, manifest = assessment_joined
    monkeypatch.setattr(current, "ROOT", root)
    monkeypatch.setattr(current.binding, "__file__", str(root / "benchmark_tools/native12_composed_review_binding.py"))
    run = dict(index=12, cell="p1_c1_r1", dataset="qfo_corrected", repeat=0, genes=3,
        inputs=stage["input_fastas"], output_root=str(root / "attempt"))
    plan_ref = put(root / "plan.json", dict(runs=[{}] * 12 + [run]))
    stage["plan"] = plan_ref
    request = dict(schema=current.binding.executor.REQUEST_SCHEMA, plan=plan_ref, index=12,
        job_id=24036, cell="p1_c1_r1", amendment=stage["amendment"],
        source=record(current.binding.executor.__file__), original_review_translated=False,
        new_sources=[record(current.binding.reviewer.__file__)])
    request_path = root / "request12.json"
    request_ref = put(request_path, request)
    monkeypatch.setattr(current.binding.executor, "REQUEST", request_path)
    monkeypatch.setattr(current.binding, "REQUEST_SHA", request_ref["sha256"])
    stage["request"] = request_ref
    stage["native_cpu_ids"] = list(range(32))
    stage["composed_binding"]["review_producer_scheduler"] = dict(JobIDRaw="25000", State="COMPLETED",
        ExitCode="0:0", NodeList="bizon", AllocCPUS="2", ReqMem="128G")
    validator = root / "benchmark_tools/validate_native_factorial_outputs.py"
    shutil.copyfile(Path(current.__file__).parent / validator.name, validator)
    outputs = dict(schema="native12_composed_output_review_v1", source=record(current.binding.reviewer.__file__),
        semantic_validator_source=record(validator), native_outputs_validated=True, accuracy_evaluated=False,
        index=12, job_id=24036, cell="p1_c1_r1", request=request_ref, plan=plan_ref,
        amendment=stage["amendment"], allocated_ready=stage["allocated_ready"], native_cpu_ids=list(range(32)),
        gene_ownership_sha256=stage["gene_ownership_sha256"], checked_files=[stage["native_input"]],
        phylogeny=dict(native_pair_rows=2))
    output_ref = put(root / "review/output.json", outputs)
    review = review_fixture()[0]
    review.update(request=request_ref, plan=plan_ref, amendment=stage["amendment"],
        source=record(current.binding.reviewer.__file__), native_cpu_ids=list(range(32)),
        resources=dict(wall_seconds=100., cpu_seconds=1000., peak_memory_bytes=1234567),
        reviews=dict(runtime={}, resources={}, environment={}, outputs_or_failure=output_ref))
    review_ref = put(root / "review/review.json", review)
    stage["terminal_review"] = review_ref
    pairs_ref = put(Path(old_ref["path"]), stage)
    context = (request, {}, {}, run, review, outputs, "native", {"fixture_native_scheduler": True}, [],
        deepcopy(stage["composed_binding"]))
    monkeypatch.setattr(runner, "native_binding", lambda *args: context)
    new_joined = (root, pairs_ref, stage, manifest)
    execution = execute(new_joined, monkeypatch)
    report = admit_fixture(new_joined, monkeypatch)
    return dict(root=root, run=run, plan_ref=plan_ref, stage=stage, request=request, review=review,
        outputs=outputs, output_ref=output_ref, report=report, execution=execution, pairs_ref=pairs_ref)


def test_real_endpoint_arithmetic_accepts_only_new_admitted_native_row(joined):
    evidence = []
    row = current.extract_final(joined["report"], joined["run"], joined["plan_ref"], evidence)
    assert row["status"] == "supplied_composed_native_admission"
    assert row["index"] == 12 and row["native_job_id"] == 24036
    assert row["accuracy_admitted"] is True and row["scientific_timings_admitted"] is False
    assert row["native_cpu_ids"] == list(range(32))
    assert row["secondary_mean"] == sum(row["scores"].values()) / 6
    assert set(row["scores"]) == set(current.original.ENDPOINTS)
    assert row["endpoint_details"]["SwissTrees"]["statistic"] == "F1"
    assert row["endpoint_details"]["GO"]["statistic"] != "F1"
    assert row["assessment_resource_limits"]["memory_bytes"] == 128 * 1024 ** 3
    assert joined["pairs_ref"] in evidence


@pytest.mark.parametrize("field,value", [("schema", "allocated_native_factorial_qfo_admission_v1"),
    ("status", "pending"), ("accuracy_admitted", False), ("native_index", 11),
    ("native_job_id", 23985), ("cell", "p1_c1_r0"), ("original_review_translated", True),
    ("next_identity_authorized", True), ("publication_ready", True), ("source", {})])
def test_other_schema_pending_or_broadened_admission_refuses(joined, field, value):
    joined["report"][field] = value
    with pytest.raises(ValueError):
        current.extract_final(joined["report"], joined["run"], joined["plan_ref"], [])


@pytest.mark.parametrize("field,value", [("phylogeny", {"native_pair_rows": 3}),
    ("schema", "allocated_native_factorial_output_review_v1"), ("accuracy_evaluated", True),
    ("native_cpu_ids", [1]), ("gene_ownership_sha256", "wrong"), ("semantic_validator_source", {})])
def test_changed_resealed_scientific_metadata_refuses(joined, field, value):
    joined["outputs"][field] = value
    output_ref = put(Path(joined["output_ref"]["path"]), joined["outputs"])
    joined["review"]["reviews"]["outputs_or_failure"] = output_ref
    review_ref = put(Path(joined["stage"]["terminal_review"]["path"]), joined["review"])
    joined["stage"]["terminal_review"] = review_ref
    joined["report"]["conversion"] = deepcopy(joined["stage"])
    joined["report"]["pairs_manifest"] = put(Path(joined["pairs_ref"]["path"]), joined["stage"])
    with pytest.raises(ValueError):
        current.extract_final(joined["report"], joined["run"], joined["plan_ref"], [])


def test_secondary_mean_cannot_be_changed_after_admission(joined):
    joined["report"]["assessment"]["secondary_six_metric_mean"] = 0.
    with pytest.raises(ValueError, match="secondary mean"):
        current.extract_final(joined["report"], joined["run"], joined["plan_ref"], [])


def baseline_fixture(tmp_path, monkeypatch, plan_ref):
    monkeypatch.setattr(current, "ROOT", tmp_path)
    cells = current.original.CELLS
    rows = []
    for index, cell in enumerate(cells, start=6):
        admitted = index in (6, 7, 8, 10)
        rows.append(dict(index=index, cell=cell, accuracy_admitted=admitted,
            status="supplied_native_admission" if admitted else "no_supplied_native_admission",
            measurement_status="fixture", scores={endpoint: .5 if admitted else None for endpoint in current.original.ENDPOINTS},
            secondary_mean=.5 if admitted else None, relation_coverage=.5 if admitted else None,
            prediction_semantics="fixture"))
    baseline = dict(schema="allocated_native_qfo_scientific_reporting_snapshot_v1", source=record(current.prior.__file__),
        supplied_admissions=4, rows=rows, plan=plan_ref, outputs=[], timing_disclosure="Shared-host fixture.")
    baseline_path = tmp_path / "baseline.json"
    baseline_ref = put(baseline_path, baseline)
    failure = dict(schema="composed_native_qfo_failure_reporting_v1", new_scientific_admission=False,
        frozen_manuscript_replaced=False, source=record(current.__file__), outputs=[], references={},
        row=dict(native_index=11, cell="p1_c1_r0", scoring_status="OUT_OF_MEMORY", accuracy_admitted=False,
            admitted_endpoint_score_count=0, secondary_six_metric_mean=None, native_job_id=23985,
            conversion_job_id=24035, assessment_job_id=24038, inference_status="native_reviewed",
            resources=dict(wall_seconds=1., cpu_seconds=2., peak_memory_bytes=123), resource_scopes=current.original.SCOPES,
            submitted_pairs=10, input_accessions=100, relation_accessions=20, relation_coverage=.2,
            prediction_semantics="cross-species group-derived clique pairs"))
    failure_path = tmp_path / "failure.json"
    failure_ref = put(failure_path, failure)
    monkeypatch.setattr(current, "BASELINE", baseline_path)
    monkeypatch.setattr(current, "BASELINE_SHA", baseline_ref["sha256"])
    monkeypatch.setattr(current, "FAILURE", failure_path)
    monkeypatch.setattr(current, "FAILURE_SHA", failure_ref["sha256"])
    return baseline, failure


def test_no_final_admission_preserves_four_rows_and_null_failed_scores(tmp_path, monkeypatch):
    plan_ref = put(tmp_path / "plan.json", dict(runs=[{}] * 13))
    baseline, failure = baseline_fixture(tmp_path, monkeypatch, plan_ref)
    result = current.collect()
    assert result["schema"] == current.SCHEMA and result["supplied_admissions"] == 4
    assert result["final_admission"] is None
    for offset in (0, 1, 2, 4):
        assert result["rows"][offset] == baseline["rows"][offset]
    failed = result["rows"][5]
    assert failed["status"] == "retained_composed_scoring_failure"
    assert failed["accuracy_admitted"] is False and failed["secondary_mean"] is None
    assert all(score is None for score in failed["scores"].values())
    assert result["rows"][6] == baseline["rows"][6]
    assert result["new_scoring_or_admission"] is False


def test_successor_collect_requires_actual_admission_not_unscored_prediction(joined, monkeypatch):
    baseline, failure = baseline_fixture(joined["root"], monkeypatch, joined["plan_ref"])
    admission_path = current.admission.DESTINATION / "results.json"
    admission_ref = record(admission_path)
    result = current.collect(admission_ref)
    assert result["supplied_admissions"] == 5
    assert result["rows"][6]["admission"] == admission_ref
    assert result["rows"][6]["status"] == "supplied_composed_native_admission"
    assert result["rows"][5]["accuracy_admitted"] is False
    for offset in (0, 1, 2, 4):
        assert result["rows"][offset] == baseline["rows"][offset]


def test_new_table_explicitly_labels_failed_scores_unavailable(tmp_path, monkeypatch):
    plan_ref = put(tmp_path / "plan.json", dict(runs=[{}] * 13))
    baseline_fixture(tmp_path, monkeypatch, plan_ref)
    output = tmp_path / "successor"
    ref = current.export(output)
    report = json.loads(Path(ref["path"]).read_text())
    assert report["new_scoring_or_admission"] is False and report["publication_ready"] is False
    table = (output / "scores.tsv").read_text()
    failed = next(line for line in table.splitlines() if "p1_c1_r0" in line)
    assert failed.count("Unavailable") == 7
    assert "GO Similarity" in table and "EC Similarity" in table
    with pytest.raises(ValueError):
        current.export(output)

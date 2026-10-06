"""Reporting recovered accuracy must preserve native measurement failure."""

import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import export_native_qfo_scientific_scores as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_export_native_qfo_factorial_scores import fixture, write
from tests.unit.test_prepare_measurement_failed_native_qfo_pairs import recovered_fixture


def recovered(tmp_path, changes=None):
    changes = changes or {}
    def change(key, value):
        if key in changes: changes[key](value)
        return value
    plan, ref = fixture(tmp_path, index=7)
    report = json.loads(Path(ref["path"]).read_text())
    stage = report["conversion"]
    participant = "ohmm_qfo_recovered_native_" + stage["cell"]
    review, _, _, _ = recovered_fixture()
    review.update(request=stage["request"], plan=plan, index=7, job_id=10)
    review_ref = write(tmp_path / "recovery.json", change("review", review))
    stage.pop("terminal_review")
    stage.update(schema="measurement_failed_native_qfo_conversion_v1",
        status="measurement_failed_native_qfo_pairs_prepared_unscored", scientific_recovery=review_ref,
        participant=participant, resources=None, original_native_scheduler_success=False,
        scientific_timings_admitted=False, eligible_for_timing_comparison=False,
        pairs=dict(bytes=123, sha256="e" * 64), filtered_pairs=dict(bytes=123, sha256="e" * 64))
    stage = change("stage", stage)
    pairs_ref = write(tmp_path / "pairs.json", stage)
    preflight = json.loads(Path(report["preflight"]["path"]).read_text())
    preflight.update(schema="measurement_failed_native_qfo_execution_v1", stage=stage, pairs_manifest=pairs_ref,
        original_native_scheduler_success=False, scientific_timings_admitted=False, eligible_for_timing_comparison=False,
        native_inference_reexecuted=False, automatic_retry=False, next_identity_authorized=False,
        publication_ready=False, resources=None, native_scheduler=dict(verified=dict(State="FAILED", ExitCode="1:0")))
    preflight = change("preflight", preflight)
    preflight_ref = write(tmp_path / "preflight.json", preflight)
    execution = dict(preflight, status="process_succeeded_pending_independent_admission", exit_code=0)
    execution_ref = write(tmp_path / "execution.json", change("execution", execution))
    report.update(schema="measurement_failed_native_qfo_admission_v1",
        status="measurement_failed_native_qfo_assessment_admitted", conversion=stage, participant=participant,
        pairs_manifest=pairs_ref, scientific_recovery=review_ref, execution_report=execution_ref, preflight=preflight_ref,
        source=record(Path(module.__file__).with_name("admit_measurement_failed_native_qfo_assessment.py")),
        original_native_scheduler_success=False, scientific_timings_admitted=False, eligible_for_timing_comparison=False,
        native_inference_reexecuted=False, resources=None, native_scheduler=dict(verified=dict(State="FAILED", ExitCode="1:0")))
    report["assessment"]["participant"] = participant
    for endpoint in report["assessment"]["endpoints"].values():
        endpoint["native_participant"]["participant_id"] = participant
    report = change("admission", report)
    return plan, write(tmp_path / "admission.json", report)


def collect(refs):
    plan, ref = refs
    return module.collect(plan["path"], plan["sha256"], [], [(ref["path"], ref["sha256"])])


def test_recovered_accuracy_keeps_timing_null_and_endpoints_separate(tmp_path):
    report = collect(recovered(tmp_path))
    row = report["rows"][1]
    assert row["scores"]["VGNC"] == pytest.approx(.6)
    assert row["endpoint_details"]["VGNC"]["precision"] == .75
    assert row["endpoint_details"]["VGNC"]["recall"] == .5
    assert row["endpoint_details"]["GO"]["assessed_relations"] == 100
    assert row["resources"] is None and row["accuracy_admitted"] is True
    assert row["timing_eligible"] is row["timing_admitted"] is False
    assert row["measurement_status"] == "failed_timing_scientific_outputs_recovered"
    assert report["supplied_admissions"] == report["supplied_recovered_admissions"] == 1
    assert report["rows"][0]["scores"]["VGNC"] is None
    assert report["new_scoring_or_admission"] is report["publication_ready"] is False


@pytest.mark.parametrize("target,key,value", [
    ("admission", "schema", "full_native_factorial_qfo_admission_v1"),
    ("admission", "status", "validating"), ("admission", "accuracy_admitted", False),
    ("admission", "publication_ready", True), ("admission", "next_identity_authorized", True),
    ("admission", "automatic_retry", True), ("admission", "native_inference_reexecuted", True),
    ("admission", "original_native_scheduler_success", True), ("admission", "scientific_timings_admitted", True),
    ("admission", "eligible_for_timing_comparison", True), ("admission", "resources", {}),
    ("admission", "native_index", True), ("admission", "native_job_id", True),
    ("admission", "cell", "p0_c0_r0"), ("admission", "scientific_recovery", {}),
    ("stage", "schema", "full_native_factorial_qfo_conversion_v1"), ("stage", "resources", {}),
    ("stage", "conversion_kind", "group"), ("stage", "removed_mapping_pairs", 1),
    ("stage", "native_job_id", 11), ("stage", "participant", "ohmm_qfo_full_native_p0_c0_r1"),
    ("review", "scheduler_state", "COMPLETED"), ("review", "scheduler_exit_code", "0:0"),
    ("review", "resources", {}), ("review", "native_outputs_validated", False),
    ("execution", "schema", "full_native_factorial_qfo_execution_v1"), ("execution", "exit_code", True),
    ("execution", "status", "running"), ("execution", "cell", "p0_c0_r0"),
    ("execution", "resources", {}), ("execution", "scientific_timings_admitted", True),
    ("preflight", "status", "unrun"),
])
def test_relabelled_or_mismatched_reports_refused(tmp_path, target, key, value):
    with pytest.raises(ValueError):
        collect(recovered(tmp_path, {target: lambda x: x.update({key: value})}))


@pytest.mark.parametrize("field,value", [("metric_x", True), ("metric_x", -.1), ("metric_x", 1.1),
    ("metric_y", float("nan")), ("metric_y", float("inf")), ("metric_y", 1.1),
    ("stderr_x", -.01), ("stderr_y", True)])
def test_invalid_endpoint_numbers_refused(tmp_path, field, value):
    with pytest.raises(ValueError):
        collect(recovered(tmp_path, {"admission": lambda x:
            x["assessment"]["endpoints"]["VGNC"]["native_participant"].update({field: value})}))


@pytest.mark.parametrize("change", ["mean", "score", "axes", "semantics", "participant", "missing", "fas", "scheduler"])
def test_endpoint_arithmetic_and_native_failure_required(tmp_path, change):
    def modify(report):
        if change == "mean": report["assessment"]["secondary_six_metric_mean"] = 0
        elif change == "missing": report["assessment"]["endpoints"].pop("EC")
        elif change == "fas": report["fas_sample"]["sample_membership_verified"] = False
        elif change == "scheduler": report["native_scheduler"]["verified"]["State"] = "COMPLETED"
        else:
            endpoint = report["assessment"]["endpoints"]["VGNC"]
            if change == "score": endpoint["score"] = 0
            elif change == "axes": endpoint["axes"]["x_axis"] = "PPV"
            elif change == "semantics": endpoint["score_semantics"] = "mean"
            else: endpoint["native_participant"]["participant_id"] = "old_cached"
    with pytest.raises(ValueError): collect(recovered(tmp_path, {"admission": modify}))


def test_mixed_normal_and_recovered_admissions(tmp_path):
    a, b = tmp_path / "normal", tmp_path / "recovered"
    a.mkdir();b.mkdir()
    plan, normal = fixture(a)
    recovered_plan, recovery = recovered(b)
    # Both immutable fixtures deliberately share one identical plan reference.
    values = json.loads(Path(recovery["path"]).read_text())
    stage = values["conversion"]
    request = json.loads(Path(stage["request"]["path"]).read_text())
    request["plan"] = plan
    request_ref = write(Path(stage["request"]["path"]), request)
    review = json.loads(Path(stage["scientific_recovery"]["path"]).read_text())
    review.update(plan=plan, request=request_ref)
    recovery_ref = write(Path(stage["scientific_recovery"]["path"]), review)
    stage.update(plan=plan, request=request_ref, scientific_recovery=recovery_ref)
    pairs = write(Path(values["pairs_manifest"]["path"]), stage)
    values.update(conversion=stage, pairs_manifest=pairs, scientific_recovery=recovery_ref)
    preflight = json.loads(Path(values["preflight"]["path"]).read_text())
    preflight.update(stage=stage, pairs_manifest=pairs)
    values["preflight"] = write(Path(values["preflight"]["path"]), preflight)
    values["execution_report"] = write(Path(values["execution_report"]["path"]),
        dict(preflight, status="process_succeeded_pending_independent_admission", exit_code=0))
    recovery = write(Path(recovery["path"]), values)
    report = module.collect(plan["path"], plan["sha256"], [(normal["path"], normal["sha256"])],
        [(recovery["path"], recovery["sha256"])])
    assert report["supplied_admissions"] == 2 and report["supplied_recovered_admissions"] == 1
    assert report["rows"][0]["status"] == "supplied_native_admission"
    assert report["rows"][1]["status"] == "supplied_recovered_scientific_admission"


def test_duplicates_and_old_exporter_reject_recovery(tmp_path):
    plan, ref = recovered(tmp_path)
    admissions = [(ref["path"], ref["sha256"])]
    with pytest.raises(ValueError): module.collect(plan["path"], plan["sha256"], [], admissions * 2)
    with pytest.raises(ValueError): module.original.collect(plan["path"], plan["sha256"], admissions)


def test_export_no_overwrite_missing_values_and_status_columns(tmp_path):
    plan, ref = recovered(tmp_path)
    output = tmp_path / "export"
    report = module.export(plan["path"], plan["sha256"], [], [(ref["path"], ref["sha256"])], output)
    with (output / "scores.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert rows[0]["VGNC F1"] == "" and rows[1]["VGNC F1"] == "0.6"
    assert rows[1]["Measurement status"] == "failed_timing_scientific_outputs_recovered"
    assert "Unavailable" in (output / "scores.md").read_text()
    assert report["recovered_inference_resources_admitted"] is False
    with pytest.raises(ValueError): module.export(plan["path"], plan["sha256"], [], [], output)


def test_same_cell_cannot_be_supplied_as_normal_and_recovered(tmp_path):
    a, b = tmp_path / "normal", tmp_path / "recovered"
    a.mkdir();b.mkdir()
    plan, normal = fixture(a, index=7)
    _, recovery = recovered(b)
    with pytest.raises(ValueError, match="duplicate"):
        module.collect(plan["path"], plan["sha256"], [(normal["path"], normal["sha256"])],
            [(recovery["path"], recovery["sha256"])])


@pytest.mark.parametrize("key,value", [("input_accessions", 9), ("accessions_in_any_pair", True),
    ("accessions_in_any_pair", 8), ("fraction_inputs_in_any_pair", .5), ("fraction_inputs_in_any_pair", True)])
def test_relation_coverage_denominator_and_pair_bound_required(tmp_path, key, value):
    with pytest.raises(ValueError):
        collect(recovered(tmp_path, {"stage": lambda x: x["pair_coverage"].update({key: value})}))


def test_changed_supplied_receipt_not_reported(tmp_path):
    plan, ref = recovered(tmp_path)
    Path(ref["path"]).write_text("{}")
    with pytest.raises(ValueError):
        module.collect(plan["path"], plan["sha256"], [], [(ref["path"], ref["sha256"])])


def test_noninteger_similarity_assessed_count_not_reported(tmp_path):
    with pytest.raises(ValueError, match="Noninteger"):
        collect(recovered(tmp_path, {"admission": lambda x:
            x["assessment"]["endpoints"]["GO"]["native_participant"].update(metric_x=1.5)}))

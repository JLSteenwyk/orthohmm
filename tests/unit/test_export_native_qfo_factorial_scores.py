import copy
import csv
import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import export_native_qfo_factorial_scores as report
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True) + "\n")
    return record(path)


def fixture(tmp_path, index=6, changes=None):
    changes = changes or {}

    def changed(name, value):
        if name in changes:
            changes[name](value)
        return value

    runs = [dict(index=i, dataset="orthobench", cell="not_used", genes=10) for i in range(6)]
    runs += [dict(index=i + 6, dataset="qfo_corrected", cell=cell, genes=10)
             for i, cell in enumerate(report.CELLS)]
    plan = write(tmp_path / "plan.json", changed("plan", dict(runs=runs)))
    cell = runs[index]["cell"]
    request = write(tmp_path / "request.json", changed("request", dict(plan=plan, index=index, job_id=10)))
    review = write(tmp_path / "review.json", changed("review", dict(
        schema="native_factorial_terminal_review_v1", request=request, plan=plan, index=index,
        cell=cell, dataset="qfo_corrected", job_id=10, status="native_success", scheduler_state="COMPLETED",
        scheduler_exit_code="0:0", execution_scope="shared_host_matched_resources", resource_scopes=report.SCOPES,
        uncontended_timing=False, terminal_reviewed=True, native_outputs_validated=True,
        primary_resources_replayed=True, shared_host_resources_reviewed=True)))
    kind = "native" if cell.endswith("r1") else "group"
    semantics = "native phylogenetically inferred pairs" if kind == "native" else "cross-species group-derived clique pairs"
    participant = "ohmm_qfo_full_native_" + cell
    stage = changed("stage", dict(schema="full_native_factorial_qfo_conversion_v1",
        status="full_native_factorial_qfo_pairs_prepared_unscored", plan=plan, request=request,
        terminal_review=review, native_index=index, native_job_id=10, job_id="20", cell=cell,
        conversion_kind=kind, semantics=semantics, participant=participant, total_pairs=3, retained_pairs=3,
        expected_pairs=3, removed_mapping_pairs=0, empty_predictions=False, accuracy_evaluated=False,
        publication_ready=False, native_inference_reexecuted=False, next_identity_authorized=False,
        automatic_retry=False, pair_coverage=dict(input_accessions=10, accessions_in_any_pair=4,
                                                 pair_rows=3, fraction_inputs_in_any_pair=.4)))
    pairs = write(tmp_path / "pairs.json", stage)
    preflight = changed("preflight", dict(schema="full_native_factorial_qfo_execution_v1", status="running",
        job_id="30", native_index=index, native_job_id=10, cell=cell, stage=stage,
        pairs_manifest=pairs, accuracy_admitted=False))
    preflight_ref = write(tmp_path / "preflight.json", preflight)
    execution = copy.deepcopy(preflight)
    execution.update(status="process_succeeded_pending_independent_admission", exit_code=0)
    execution_ref = write(tmp_path / "execution.json", changed("execution", execution))
    endpoints = {}
    for name in report.ENDPOINTS:
        is_f1 = name in report.F1_ENDPOINTS
        x, y = (.5, .75) if is_f1 else (100, {"GO": .2, "EC": .4, "FAS": .8}[name])
        endpoints[name] = dict(score=report.harmonic_mean(x, y) if is_f1 else y,
            score_semantics="harmonic mean of native TPR and PPV" if is_f1 else report.AXES[name][1],
            axes=dict(x_axis=report.AXES[name][0], y_axis=report.AXES[name][1]),
            native_participant=dict(participant_id=participant, metric_x=x, metric_y=y, stderr_x=.01, stderr_y=.02))
    admitted = changed("admission", dict(schema="full_native_factorial_qfo_admission_v1",
        status="full_native_factorial_qfo_assessment_admitted", accuracy_admitted=True,
        publication_ready=False, next_identity_authorized=False, automatic_retry=False,
        native_index=index, native_job_id=10, cell=cell, participant=participant,
        pairs_manifest=pairs, conversion=stage, execution_report=execution_ref, preflight=preflight_ref,
        scheduler=dict(JobIDRaw="30", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="8"),
        conversion_scheduler=dict(JobIDRaw="20", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2"),
        source=record(Path(report.__file__).with_name("admit_native_factorial_qfo_assessment.py")),
        assessment=dict(participant=participant, endpoints=endpoints,
                        secondary_six_metric_mean=sum(e["score"] for e in endpoints.values()) / 6),
        fas_sample=dict(sample_membership_verified=True), fas_protocol=dict(population="eligible", sampling="unseeded")))
    admission = write(tmp_path / "admission.json", admitted)
    return plan, admission


def collect(refs):
    plan, admission = refs
    return report.collect(plan["path"], plan["sha256"], [(admission["path"], admission["sha256"])])


def test_admitted_scores_keep_endpoint_semantics_and_coverage(tmp_path):
    result = collect(fixture(tmp_path))
    assert len(result["rows"]) == 7 and result["supplied_admissions"] == 1
    row = result["rows"][0]
    assert row["scores"]["VGNC"] == pytest.approx(.6)
    assert row["endpoint_details"]["VGNC"]["precision"] == .75
    assert row["endpoint_details"]["VGNC"]["recall"] == .5
    assert row["endpoint_details"]["GO"]["statistic"] == "avg Schlicker"
    assert row["endpoint_details"]["GO"]["assessed_relations"] == 100
    assert row["submitted_pairs"] == 3 and row["input_accessions"] == 10
    assert row["relation_accessions"] == 4 and row["relation_coverage"] == .4
    assert row["secondary_mean"] == pytest.approx(3.2 / 6)
    assert result["rows"][1]["scores"]["VGNC"] is None
    assert result["new_scoring_or_admission"] is False and result["publication_ready"] is False


def test_resolved_predictions_not_reported_as_group_cliques(tmp_path):
    row = collect(fixture(tmp_path, index=7))["rows"][1]
    assert row["prediction_semantics"] == "native phylogenetically inferred pairs"


@pytest.mark.parametrize("target,key,value", [
    ("admission", "status", "validating"), ("admission", "schema", "corrected_factorial_assessment_admitted"),
    ("admission", "accuracy_admitted", False), ("admission", "publication_ready", True),
    ("admission", "automatic_retry", True), ("admission", "next_identity_authorized", True),
    ("admission", "native_index", True), ("admission", "native_index", 5),
    ("admission", "native_job_id", True), ("admission", "native_job_id", 11),
    ("admission", "cell", "p1_c0_r0"),
    ("stage", "conversion_kind", "native"), ("stage", "semantics", "pre-clustering edges"),
    ("stage", "total_pairs", 4), ("stage", "removed_mapping_pairs", 1),
    ("stage", "empty_predictions", True), ("stage", "accuracy_evaluated", True),
    ("stage", "native_inference_reexecuted", True), ("stage", "status", "preparing_unscored"),
    ("request", "job_id", 11), ("review", "job_id", 11), ("review", "uncontended_timing", True),
    ("review", "status", "native_failure_retained"), ("review", "scheduler_state", "RUNNING"),
    ("review", "native_outputs_validated", False), ("review", "primary_resources_replayed", False),
    ("execution", "exit_code", 1), ("execution", "status", "running"),
    ("execution", "cell", "p0_c1_r0"), ("execution", "accuracy_admitted", True),
    ("preflight", "status", "not_started"),
])
def test_invalid_admission_or_binding_refused(tmp_path, target, key, value):
    refs = fixture(tmp_path, changes={target: lambda v: v.update({key: value})})
    with pytest.raises(ValueError):
        collect(refs)


@pytest.mark.parametrize("field,value", [
    ("metric_x", True), ("metric_x", -.1), ("metric_x", 1.1),
    ("metric_y", float("nan")), ("metric_y", float("inf")),
    ("metric_y", 1.1), ("stderr_x", -.01), ("stderr_y", True),
])
def test_invalid_endpoint_number_refused(tmp_path, field, value):
    refs = fixture(tmp_path, changes={"admission": lambda v:
        v["assessment"]["endpoints"]["VGNC"]["native_participant"].update({field: value})})
    with pytest.raises(ValueError):
        collect(refs)


@pytest.mark.parametrize("target,key,value", [
    ("pair_coverage", "input_accessions", 9), ("pair_coverage", "accessions_in_any_pair", 11),
    ("pair_coverage", "accessions_in_any_pair", True),
    ("pair_coverage", "fraction_inputs_in_any_pair", .5), ("pair_coverage", "pair_rows", 4),
    ("pair_coverage", "fraction_inputs_in_any_pair", True), ("pair_coverage", "pair_rows", True),
])
def test_wrong_relation_coverage_refused(tmp_path, target, key, value):
    refs = fixture(tmp_path, changes={"stage": lambda v: v[target].update({key: value})})
    with pytest.raises(ValueError, match="coverage"):
        collect(refs)


def test_no_pairs_cannot_have_relation_coverage(tmp_path):
    def change(stage):
        stage.update(total_pairs=0, retained_pairs=0, expected_pairs=0, empty_predictions=True)
        stage["pair_coverage"].update(pair_rows=0)
    with pytest.raises(ValueError, match="coverage"):
        collect(fixture(tmp_path, changes={"stage": change}))


def test_admitted_zero_endpoint_is_not_replaced_with_missing(tmp_path):
    def change(admission):
        endpoint = admission["assessment"]["endpoints"]["VGNC"]
        endpoint["native_participant"].update(metric_x=0., metric_y=0.)
        endpoint["score"] = 0.
        admission["assessment"]["secondary_six_metric_mean"] = sum(
            v["score"] for v in admission["assessment"]["endpoints"].values()) / 6
    result = collect(fixture(tmp_path, changes={"admission": change}))
    assert result["rows"][0]["scores"]["VGNC"] == 0.
    assert result["rows"][1]["scores"]["VGNC"] is None


def test_boolean_plan_index_refused(tmp_path):
    with pytest.raises(ValueError, match="plan"):
        collect(fixture(tmp_path, changes={"plan": lambda p: p["runs"][1].update(index=True)}))


@pytest.mark.parametrize("change", [
    lambda v: v["assessment"]["endpoints"].pop("GO"),
    lambda v: v["assessment"]["endpoints"]["GO"].update(score=.3),
    lambda v: v["assessment"].update(secondary_six_metric_mean=.1),
    lambda v: v["assessment"]["endpoints"]["GO"]["native_participant"].update(metric_x=100.5),
    lambda v: v["assessment"]["endpoints"]["GO"].update(score_semantics="F1"),
    lambda v: v["assessment"]["endpoints"]["VGNC"]["axes"].update(x_axis="PPV"),
    lambda v: v["assessment"]["endpoints"]["GO"]["native_participant"].update(participant_id="cached"),
    lambda v: v["scheduler"].update(State="RUNNING"),
    lambda v: v["scheduler"].update(AllocCPUS="2"),
    lambda v: v["conversion_scheduler"].update(ExitCode="1:0"),
    lambda v: v["source"].update(sha256="0" * 64),
    lambda v: v["conversion"].update(cell="p0_c1_r1"),
    lambda v: v["fas_sample"].update(sample_membership_verified=False),
])
def test_arithmetic_identity_and_scheduler_errors_refused(tmp_path, change):
    refs = fixture(tmp_path, changes={"admission": change})
    with pytest.raises(ValueError):
        collect(refs)


def test_checksums_duplicates_and_symlinks_refused(tmp_path):
    plan, admission = fixture(tmp_path)
    arg = (admission["path"], admission["sha256"])
    with pytest.raises(ValueError, match="duplicate"):
        report.collect(plan["path"], plan["sha256"], [arg, arg])
    with pytest.raises(ValueError, match="checksum"):
        report.collect(plan["path"], "0" * 64, [])
    link = tmp_path / "link.json"
    link.symlink_to(admission["path"])
    with pytest.raises(ValueError, match="Nonregular"):
        report.collect(plan["path"], plan["sha256"], [(str(link), admission["sha256"])])
    Path(admission["path"]).write_text("{}\n")
    with pytest.raises(ValueError, match="checksum"):
        report.collect(plan["path"], plan["sha256"], [arg])


def test_empty_export_keeps_unavailable_instead_of_zero(tmp_path):
    plan, _ = fixture(tmp_path)
    result = report.export(plan["path"], plan["sha256"], [], tmp_path / "empty")
    assert result["supplied_admissions"] == 0
    assert all(r["relation_coverage"] is None for r in result["rows"])
    assert "Unavailable" in (tmp_path / "empty/scores.md").read_text()


def test_export_writes_exact_tsv_labels_and_refuses_overwrite(tmp_path):
    plan, admission = fixture(tmp_path)
    output = tmp_path / "table"
    result = report.export(plan["path"], plan["sha256"], [(admission["path"], admission["sha256"])], output)
    with (output / "scores.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 7 and rows[0]["VGNC F1"] == "0.6" and rows[0]["Relation coverage"] == "0.4"
    assert rows[1]["VGNC F1"] == "" and rows[1]["Secondary mean"] == ""
    text = (output / "scores.md").read_text()
    assert "GO similarity" in text and "Secondary mean" in text and report.DISCLOSURE in text
    assert "not official QfO F1" in text and "P1C0R0" in text
    assert json.loads((output / "report.json").read_text())["rows"] == result["rows"]
    with pytest.raises(ValueError, match="already exists"):
        report.export(plan["path"], plan["sha256"], [], output)


def test_failed_binding_creates_no_report_directory(tmp_path):
    plan, admission = fixture(tmp_path, changes={"admission": lambda v: v.update(accuracy_admitted=False)})
    output = tmp_path / "table"
    with pytest.raises(ValueError):
        report.export(plan["path"], plan["sha256"], [(admission["path"], admission["sha256"])], output)
    assert not output.exists()


def test_cli_uses_supplied_native_admission_only(tmp_path):
    plan, admission = fixture(tmp_path)
    output = tmp_path / "cli"
    p = subprocess.run([sys.executable, "-B", report.__file__, "--plan", plan["path"],
                        "--plan-sha256", plan["sha256"], "--admission", admission["path"],
                        admission["sha256"], "--output", str(output)], capture_output=True, text=True, check=True)
    assert json.loads(p.stdout) == dict(supplied_admissions=1, rows=7)

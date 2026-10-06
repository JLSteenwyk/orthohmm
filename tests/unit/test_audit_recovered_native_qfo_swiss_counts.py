import copy
import gzip
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import audit_native_qfo_swiss_counts as ordinary
from benchmark_tools import audit_qfo_factorial_swiss as retained_audit
from benchmark_tools import audit_recovered_native_qfo_swiss_counts as module
from benchmark_tools import bootstrap_qfo_factorial as bootstrap
from benchmark_tools import export_native_qfo_scientific_scores as reporter
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_audit_qfo_factorial_swiss import synthetic
from tests.unit.test_export_native_qfo_factorial_scores import fixture as normal_fixture
from tests.unit.test_export_native_qfo_scientific_scores import recovered, write


def fixture(tmp_path, monkeypatch, problem=None, different_counts=False):
    raw_root = tmp_path / "retained_raw"
    raw_root.mkdir()
    entries, baseline = synthetic(raw_root)
    reference = tmp_path / "reference"
    reference.write_text("synthetic reference\n")
    baseline["reference"] = record(reference)
    monkeypatch.setattr(retained_audit, "REFERENCE_SHA", baseline["reference"]["sha256"])
    monkeypatch.setattr(bootstrap, "REFERENCE_SHA", baseline["reference"]["sha256"])
    retained = retained_audit.assemble(entries, baseline)
    retained.update(status="corrected_qfo_factorial_swiss_counts_verified", uncertainty_admitted=False)
    retained_ref = write(tmp_path / "retained.json", retained)
    monkeypatch.setattr(ordinary, "RETAINED_COUNTS_SHA", retained_ref["sha256"])
    native = tmp_path / "native"
    native.mkdir()
    plan, ref = recovered(native)
    report = json.loads(Path(ref["path"]).read_text())
    entry = copy.deepcopy(entries[2 if different_counts else 1])
    raw = native / "SwissTrees/native.raw.txt.gz"
    raw.parent.mkdir()
    raw.write_bytes(Path(entry["raw_file"]["path"]).read_bytes())
    if problem in ("truth", "members", "duplicate"):
        with gzip.open(raw, "rt") as stream:
            lines = stream.readlines()
        if problem == "truth":
            a, b = lines[1].rstrip().split("\t"), lines[2].rstrip().split("\t")
            a[-1], b[-1] = b[-1], a[-1]
            lines[1], lines[2] = "\t".join(a) + "\n", "\t".join(b) + "\n"
        elif problem == "members":
            lines = [line.replace("family0_gene0", "family0_newgene") for line in lines]
        else:
            lines.append(lines[1])
        with gzip.open(raw, "wt") as stream:
            stream.writelines(lines)
    raw_ref = record(raw)
    assessment = entry["assessment"]
    assessment["participant"] = report["participant"]
    for metric in assessment["native_assessments"]:
        metric["participant_id"] = report["participant"]
    report["assessment"].update(assessment)
    means = {m["metrics"]["metric_id"]: m["metrics"]["value"] for m in assessment["native_assessments"]
             if m["challenge_id"] == "SwissTrees"}
    endpoint = report["assessment"]["endpoints"]["SwissTrees"]
    endpoint["native_participant"].update(metric_x=means["TPR"], metric_y=means["PPV"])
    endpoint["score"] = reporter.original.harmonic_mean(means["TPR"], means["PPV"])
    report["assessment"]["secondary_six_metric_mean"] = sum(
        e["score"] for e in report["assessment"]["endpoints"].values()) / 6
    if problem == "inventory":
        assessment["swiss_reference_families"].reverse()
    elif problem == "family_metric":
        assessment["native_assessments"][-1]["metrics"]["value"] = .99
    elif problem == "missing_metric":
        assessment["native_assessments"].pop()
    preflight = json.loads(Path(report["preflight"]["path"]).read_text())
    preflight["outputs"] = [] if problem == "missing_raw" else [raw_ref]
    if problem == "ambiguous_raw":
        preflight["outputs"] *= 2
    report["preflight"] = write(Path(report["preflight"]["path"]), preflight)
    report["execution_report"] = write(Path(report["execution_report"]["path"]),
        dict(preflight, status="process_succeeded_pending_independent_admission", exit_code=0))
    report["checked_records"] = [report["execution_report"], raw_ref]
    if problem == "unchecked_execution":
        report["checked_records"].remove(report["execution_report"])
    elif problem == "unchecked_raw":
        report["checked_records"].remove(raw_ref)
    ref = write(Path(ref["path"]), report)
    output = tmp_path / "snapshot"
    reporter.export(plan["path"], plan["sha256"], [], [(ref["path"], ref["sha256"])], output)
    return record(output / "report.json"), retained_ref, raw_ref


def run(refs):
    snapshot, retained, _ = refs
    return module.audit(snapshot["path"], snapshot["sha256"], retained["path"])


def test_admitted_recovered_counts_keep_timing_failed_and_no_intervals(tmp_path, monkeypatch):
    result = run(fixture(tmp_path, monkeypatch))
    assert len(result["families"]) == 18 and len(result["cells"]) == 1
    assert len(result["unavailable_cells"]) == 6
    assert result["successful_native_cells_not_recounted"] == []
    assert result["reference_relation_count"] == 270
    cell = result["cells"][0]
    assert cell["cell"] == "p0_c0_r1" and cell["index"] == 7
    assert cell["retained_family_records_identical"] is cell["retained_aggregate_identical"] is True
    assert cell["resources"] is None and cell["timing_admitted"] is cell["timing_eligible"] is False
    assert cell["measurement_status"] == "failed_timing_scientific_outputs_recovered"
    assert result["new_bootstrap_draws"] == 0
    for flag in ("historical_intervals_attached", "independent_confirmation", "publication_ready",
                 "new_accuracy_or_resource_admission"):
        assert result[flag] is False


def test_differing_actual_counts_not_replaced_by_cached_counts(tmp_path, monkeypatch):
    result = run(fixture(tmp_path, monkeypatch, different_counts=True))
    assert result["cells"][0]["retained_family_records_identical"] is False
    assert result["cells"][0]["differing_families"] == result["families"]


def test_mixed_snapshot_does_not_recount_successful_native_cell(tmp_path, monkeypatch):
    refs = fixture(tmp_path, monkeypatch)
    snapshot = json.loads(Path(refs[0]["path"]).read_text())
    normal_root = tmp_path / "normal"
    normal_root.mkdir()
    _, normal_ref = normal_fixture(normal_root)
    normal = json.loads(Path(normal_ref["path"]).read_text())
    stage = normal["conversion"]
    request = json.loads(Path(stage["request"]["path"]).read_text())
    request["plan"] = snapshot["plan"]
    request_ref = write(Path(stage["request"]["path"]), request)
    review = json.loads(Path(stage["terminal_review"]["path"]).read_text())
    review.update(plan=snapshot["plan"], request=request_ref)
    review_ref = write(Path(stage["terminal_review"]["path"]), review)
    stage.update(plan=snapshot["plan"], request=request_ref, terminal_review=review_ref)
    pairs_ref = write(Path(normal["pairs_manifest"]["path"]), stage)
    normal.update(conversion=stage, pairs_manifest=pairs_ref)
    preflight = json.loads(Path(normal["preflight"]["path"]).read_text())
    preflight.update(stage=stage, pairs_manifest=pairs_ref, outputs=[{"not_recounted": True}])
    normal["preflight"] = write(Path(normal["preflight"]["path"]), preflight)
    normal["execution_report"] = write(Path(normal["execution_report"]["path"]),
        dict(preflight, status="process_succeeded_pending_independent_admission", exit_code=0))
    normal_ref = write(Path(normal_ref["path"]), normal)
    recovered_ref = snapshot["rows"][1]["admission"]
    output = tmp_path / "mixed"
    reporter.export(snapshot["plan"]["path"], snapshot["plan"]["sha256"],
        [(normal_ref["path"], normal_ref["sha256"])],
        [(recovered_ref["path"], recovered_ref["sha256"])], output)
    result = run((record(output / "report.json"), refs[1], refs[2]))
    assert len(result["cells"]) == 1 and len(result["unavailable_cells"]) == 5
    assert result["successful_native_cells_not_recounted"] == ["p0_c0_r0"]


@pytest.mark.parametrize("problem", ["truth", "members", "duplicate", "inventory", "family_metric",
    "missing_metric", "missing_raw", "ambiguous_raw", "unchecked_execution", "unchecked_raw"])
def test_bad_or_unchecked_raw_rejected(tmp_path, monkeypatch, problem):
    refs = fixture(tmp_path, monkeypatch, problem=problem)
    with pytest.raises(ValueError):
        run(refs)


@pytest.mark.parametrize("field,value", [("schema", "native_qfo_reporting_snapshot_v1"), ("source", {}),
    ("publication_ready", True), ("supplied_recovered_admissions", 0), ("supplied_admissions", 7)])
def test_snapshot_schema_source_and_claims_replay(tmp_path, monkeypatch, field, value):
    refs = fixture(tmp_path, monkeypatch)
    snapshot = json.loads(Path(refs[0]["path"]).read_text())
    snapshot[field] = value
    changed = write(Path(refs[0]["path"]), snapshot)
    with pytest.raises(ValueError):
        run((changed, refs[1], refs[2]))


@pytest.mark.parametrize("target", ["raw", "retained", "snapshot_output"])
def test_post_snapshot_mutation_rejected(tmp_path, monkeypatch, target):
    refs = fixture(tmp_path, monkeypatch)
    path = Path(refs[{"raw": 2, "retained": 1}.get(target, 0)]["path"])
    if target == "snapshot_output":
        path = path.with_name("scores.tsv")
    path.write_bytes(b"changed")
    with pytest.raises(ValueError):
        run(refs)


def test_empty_snapshot_does_not_supply_cached_counts(tmp_path, monkeypatch):
    refs = fixture(tmp_path, monkeypatch)
    snapshot = json.loads(Path(refs[0]["path"]).read_text())
    output = tmp_path / "empty"
    reporter.export(snapshot["plan"]["path"], snapshot["plan"]["sha256"], [], [], output)
    with pytest.raises(ValueError, match="No supplied recovered accuracy"):
        run((record(output / "report.json"), refs[1], refs[2]))


def test_no_overwrite_before_audit(tmp_path, monkeypatch):
    output = tmp_path / "existing"
    output.write_text("retain me")
    monkeypatch.setattr(sys, "argv", ["audit", "--snapshot", "absent", "--snapshot-sha256", "unused",
        "--retained-counts", "absent", "--output", str(output)])
    with pytest.raises(ValueError, match="Output already exists"):
        module.main()
    assert output.read_text() == "retain me"

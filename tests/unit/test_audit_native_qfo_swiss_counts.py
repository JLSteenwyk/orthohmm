import copy
import gzip
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import audit_native_qfo_swiss_counts as audit
from benchmark_tools import audit_qfo_factorial_swiss as retained_audit
from benchmark_tools import bootstrap_qfo_factorial as bootstrap
from benchmark_tools import export_native_qfo_factorial_scores as reporter
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_audit_qfo_factorial_swiss import synthetic
from tests.unit.test_export_native_qfo_factorial_scores import fixture as native_fixture, write


def fixture(tmp_path, monkeypatch, problem=None, different_counts=False):
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    entries, baseline = synthetic(raw_dir)
    reference = tmp_path / "reference"
    reference.write_text("synthetic reference\n")
    baseline["reference"] = record(reference)
    monkeypatch.setattr(retained_audit, "REFERENCE_SHA", baseline["reference"]["sha256"])
    monkeypatch.setattr(bootstrap, "REFERENCE_SHA", baseline["reference"]["sha256"])
    retained = retained_audit.assemble(entries, baseline)
    retained.update(status="corrected_qfo_factorial_swiss_counts_verified", uncertainty_admitted=False)
    retained_ref = write(tmp_path / "retained.json", retained)
    monkeypatch.setattr(audit, "RETAINED_COUNTS_SHA", retained_ref["sha256"])
    native_dir = tmp_path / "native"
    native_dir.mkdir()
    plan, admission_ref = native_fixture(native_dir)
    admission = json.loads(Path(admission_ref["path"]).read_text())
    entry = copy.deepcopy(entries[2 if different_counts else 0])
    raw = native_dir / "SwissTrees/native.raw.txt.gz"
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
    native = entry["assessment"]
    native["participant"] = admission["participant"]
    for metric in native["native_assessments"]:
        metric["participant_id"] = admission["participant"]
    admission["assessment"].update(native)
    means = {m["metrics"]["metric_id"]: m["metrics"]["value"] for m in native["native_assessments"]
             if m["challenge_id"] == "SwissTrees"}
    endpoint = admission["assessment"]["endpoints"]["SwissTrees"]
    endpoint["native_participant"].update(metric_x=means["TPR"], metric_y=means["PPV"])
    endpoint["score"] = reporter.harmonic_mean(means["TPR"], means["PPV"])
    admission["assessment"]["secondary_six_metric_mean"] = sum(
        e["score"] for e in admission["assessment"]["endpoints"].values()) / 6
    if problem == "inventory":
        native["swiss_reference_families"].reverse()
    elif problem == "family_metric":
        native["native_assessments"][-1]["metrics"]["value"] = .99
    elif problem == "missing_metric":
        native["native_assessments"].pop()
    execution = json.loads(Path(admission["execution_report"]["path"]).read_text())
    execution["outputs"] = [raw_ref]
    if problem == "missing_raw":
        execution["outputs"] = []
    elif problem == "ambiguous_raw":
        execution["outputs"] *= 2
    execution_ref = write(Path(admission["execution_report"]["path"]), execution)
    admission["execution_report"] = execution_ref
    admission["checked_records"] = [execution_ref, raw_ref]
    if problem == "unchecked_execution":
        admission["checked_records"].remove(execution_ref)
    elif problem == "unchecked_raw":
        admission["checked_records"].remove(raw_ref)
    admission_ref = write(Path(admission_ref["path"]), admission)
    snapshot_dir = tmp_path / "snapshot"
    reporter.export(plan["path"], plan["sha256"], [(admission_ref["path"], admission_ref["sha256"])], snapshot_dir)
    snapshot = record(snapshot_dir / "report.json")
    return snapshot, retained_ref, raw_ref


def run(refs):
    snapshot, retained, _ = refs
    return audit.audit(snapshot["path"], snapshot["sha256"], retained["path"])


def test_actual_raw_counts_bound_to_native_admission_without_intervals(tmp_path, monkeypatch):
    result = run(fixture(tmp_path, monkeypatch))
    assert len(result["cells"]) == 1 and len(result["families"]) == 18
    assert result["reference_relation_count"] == 270
    assert len(result["unavailable_cells"]) == 6
    cell = result["cells"][0]
    assert cell["cell"] == "p0_c0_r0" and cell["index"] == 6
    assert cell["retained_family_records_identical"] is True
    assert cell["retained_aggregate_identical"] is True
    assert cell["differing_families"] == []
    assert cell["families"][0]["counts_without_prior"] == {"TP": 3, "FN": 5, "FP": 2, "TN": 5}
    assert result["new_bootstrap_draws"] == 0
    for flag in ("historical_intervals_attached", "independent_confirmation", "publication_ready",
                 "new_accuracy_or_resource_admission"):
        assert result[flag] is False


def test_changed_native_counts_are_retained_not_replaced_or_rejected(tmp_path, monkeypatch):
    result = run(fixture(tmp_path, monkeypatch, different_counts=True))
    cell = result["cells"][0]
    assert cell["retained_family_records_identical"] is False
    assert cell["differing_families"] == result["families"]
    assert cell["families"][0]["counts_without_prior"] != {"TP": 3, "FN": 5, "FP": 2, "TN": 5}
    assert result["historical_intervals_attached"] is False


@pytest.mark.parametrize("problem", ["truth", "members", "duplicate", "inventory", "family_metric",
                                    "missing_metric", "missing_raw", "ambiguous_raw",
                                    "unchecked_execution", "unchecked_raw"])
def test_changed_or_unchecked_raw_evidence_refused(tmp_path, monkeypatch, problem):
    refs = fixture(tmp_path, monkeypatch, problem=problem)
    with pytest.raises(ValueError):
        run(refs)


@pytest.mark.parametrize("target", ["raw", "retained", "snapshot_output"])
def test_post_export_tampering_refused(tmp_path, monkeypatch, target):
    refs = fixture(tmp_path, monkeypatch)
    path = Path(refs[{"raw": 2, "retained": 1}.get(target, 0)]["path"])
    if target == "snapshot_output":
        path = path.with_name("scores.tsv")
    path.write_bytes(b"changed")
    with pytest.raises(ValueError):
        run(refs)


@pytest.mark.parametrize("field,value", [("publication_ready", True), ("supplied_admissions", 7),
                                        ("schema", "cached_snapshot"), ("source", {})])
def test_snapshot_claims_must_replay(tmp_path, monkeypatch, field, value):
    refs = fixture(tmp_path, monkeypatch)
    snapshot = json.loads(Path(refs[0]["path"]).read_text())
    snapshot[field] = value
    new_ref = write(Path(refs[0]["path"]), snapshot)
    with pytest.raises(ValueError):
        run((new_ref, refs[1], refs[2]))


def test_no_native_admission_not_replaced_with_cached_counts(tmp_path, monkeypatch):
    refs = fixture(tmp_path, monkeypatch)
    original = json.loads(Path(refs[0]["path"]).read_text())
    destination = tmp_path / "empty"
    reporter.export(original["plan"]["path"], original["plan"]["sha256"], [], destination)
    with pytest.raises(ValueError, match="No supplied native admission"):
        run((record(destination / "report.json"), refs[1], refs[2]))


def test_main_refuses_existing_output_before_audit(tmp_path, monkeypatch):
    path = tmp_path / "existing.json"
    path.write_text("retain me")
    monkeypatch.setattr(sys, "argv", ["audit", "--snapshot", "absent", "--snapshot-sha256", "unused",
                                    "--retained-counts", "absent", "--output", str(path)])
    with pytest.raises(ValueError, match="Output already exists"):
        audit.main()
    assert path.read_text() == "retain me"

from pathlib import Path
from copy import deepcopy
import json

import pytest

from benchmark_tools import audit_qfo_corrected_swiss as module
from benchmark_tools.audit_qfo_swiss_counts import read_raw
from tests.unit.test_audit_qfo_factorial_swiss import synthetic


def evidence(tmp_path):
    entries, baseline = synthetic(tmp_path)
    families = baseline["families"]
    anchor = read_raw(Path(entries[0]["raw_file"]["path"]), families)
    raw = read_raw(Path(entries[2]["raw_file"]["path"]), families)
    return entries[2]["assessment"], raw, baseline, anchor


def test_fresh_counts_not_historical_counts(tmp_path):
    assessment, raw, baseline, anchor = evidence(tmp_path)
    result = module.compare(assessment, raw, baseline, anchor)
    assert raw[0] != anchor[0]
    assert result["reference_relation_count"] == 270
    assert len(result["families"]) == 18
    assert result["families"][0]["counts_without_prior"] == dict(raw[0]["family0"])
    p, r = (result["aggregate"][k] for k in ("PPV", "TPR"))
    assert result["aggregate"]["F1"] == pytest.approx(2 * p * r / (p + r))


@pytest.mark.parametrize("mutation", ["family", "reference", "status", "truth", "members", "coverage", "score", "participant"])
def test_bad_evidence_rejected(tmp_path, mutation):
    assessment, raw, baseline, anchor = evidence(tmp_path)
    if mutation == "family":
        assessment["swiss_reference_families"] = list(reversed(baseline["families"]))
    elif mutation == "reference":
        baseline["reference"]["sha256"] = "0" * 64
    elif mutation == "status":
        baseline["status"] = "unverified"
    elif mutation == "truth":
        key = next(iter(raw[1]))
        raw[1][key] = not raw[1][key]
    elif mutation == "members":
        raw[2]["family0"].add("new")
    elif mutation == "coverage":
        raw[0]["family0"]["TN"] += 1
    elif mutation == "score":
        assessment["native_assessments"][-1]["metrics"]["value"] = .999
    else:
        assessment["participant"] = "wrong"
    with pytest.raises(ValueError):
        module.compare(assessment, raw, baseline, anchor)


def audit_fixture(tmp_path, monkeypatch, recovered):
    entries, baseline = synthetic(tmp_path)

    def save(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
        return module.record(path)

    reference = tmp_path / "reference"
    reference.write_text("fixture reference")
    baseline["reference"] = module.record(reference)
    monkeypatch.setattr(module, "REFERENCE_SHA", baseline["reference"]["sha256"])
    baseline_record = save(tmp_path / "baseline.json", baseline)
    monkeypatch.setattr(module, "BASE_COUNTS_SHA", baseline_record["sha256"])
    conversion = {"participant": entries[2]["assessment"]["participant"], "fixture": True}
    pairs = save(tmp_path / "pairs.json", conversion)
    raw_path = tmp_path / "SwissTrees" / "raw.txt.gz"
    raw_path.parent.mkdir()
    raw_path.write_bytes(Path(entries[2]["raw_file"]["path"]).read_bytes())
    execution = dict(status="process_succeeded_pending_independent_admission", exit_code=0,
                     method="orthomcl", job_id="123", pairs_manifest=pairs,
                     stage=conversion, outputs=[module.record(raw_path)])
    directory = "qfo_blast_recovery_assessment_v1" if recovered else "qfo_corrected_assessment_v1/orthomcl"
    executed = save(tmp_path / directory / "results.json", execution)
    report = dict(status="recovered_orthomcl_assessment_admitted" if recovered else "admitted",
                  method="orthomcl", pairs_manifest=pairs, execution_report=executed,
                  checked_records=[executed], assessment=entries[2]["assessment"],
                  scheduler=dict(JobIDRaw="123", State="COMPLETED", ExitCode="0:0"),
                  query_coverage={"failed_queries": ["failed_fixture"]},
                  group_coverage={"fixture": True}, group_audit={"fixture": True},
                  pair_semantics="cross_species_final_group_cliques")
    # Exporter admission-contract tests cover extract; this fixture isolates raw-count plumbing.
    monkeypatch.setattr(module, "extract", lambda report, conversion: {"key": "orthomcl_1_4"})
    return report, baseline_record, save


@pytest.mark.parametrize("recovered", [False, True])
def test_audit_execution_inventory_and_recovery_provenance(tmp_path, monkeypatch, recovered):
    report, baseline, save = audit_fixture(tmp_path, monkeypatch, recovered)
    admission = save(tmp_path / "admission.json", report)
    result = module.audit(Path(admission["path"]), admission["sha256"], Path(baseline["path"]))
    assert result["reference_relation_count"] == 270
    assert len(result["families"]) == 18
    assert result["publication_ready"] is False
    assert result["uncertainty_admitted"] is False
    if recovered:
        assert result["search_recovery"] is True
        for key in ("query_coverage", "group_coverage", "group_audit", "pair_semantics"):
            assert result[key] == report[key]
    else:
        assert "search_recovery" not in result


@pytest.mark.parametrize("mutation", ["pin", "missing", "duplicate", "stage", "job", "raw_hash"])
def test_recovered_execution_mismatch_rejected(tmp_path, monkeypatch, mutation):
    report, baseline, save = audit_fixture(tmp_path, monkeypatch, True)
    if mutation == "pin":
        report["execution_report"] = deepcopy(report["execution_report"])
        report["execution_report"]["sha256"] = "0" * 64
    elif mutation == "missing":
        report["checked_records"] = []
    elif mutation == "duplicate":
        report["checked_records"] *= 2
    else:
        path = Path(report["execution_report"]["path"])
        execution = json.loads(path.read_text())
        if mutation == "stage":
            execution["stage"]["participant"] = "wrong"
        elif mutation == "job":
            execution["job_id"] = "124"
        else:
            execution["outputs"][0]["sha256"] = "0" * 64
        updated = save(path, execution)
        report["execution_report"] = updated
        report["checked_records"] = [updated]
    admission = save(tmp_path / "admission.json", report)
    with pytest.raises(ValueError):
        module.audit(Path(admission["path"]), admission["sha256"], Path(baseline["path"]))

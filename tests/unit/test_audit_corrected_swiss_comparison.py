from copy import deepcopy
from pathlib import Path

import pytest

from benchmark_tools import audit_corrected_swiss_comparison as module
from tests.unit.test_bootstrap_corrected_swiss_comparators import fixture


@pytest.mark.parametrize("recovered", [False, True])
def test_assembly_retains_recovery_without_changing_counts(monkeypatch, recovered):
    counts = fixture()
    target = counts["methods"][-1]
    comparison = dict(status="corrected_qfo_publication_comparison", publication_ready=False,
                      admitted_methods=1, replay_admission={"path": "replay", "sha256": module.REPLAY_SHA},
                      checked_records=[], methods=[
                          {"key": row["method"], "status": "not_admitted", "scores": {"SwissTrees": None}}
                          for row in counts["methods"]])
    comparison["methods"][-1].update(status="admitted", admission={"path": "admission", "sha256": "a"},
                                    conversion={"path": "conversion", "sha256": "c"},
                                    prediction_semantics="final MCL group-derived pairs",
                                    scores={"SwissTrees": target["aggregate"]["F1"]})
    baseline = {key: counts[key] for key in ("families", "reference", "shared_represented_genes")}
    baseline["stages"] = [{"raw_file": {"path": "anchor"}}]
    report = dict(status="recovered_orthomcl_assessment_admitted" if recovered else "original_admission",
                  pairs_manifest=comparison["methods"][-1]["conversion"])
    audited = dict(target, method_key=target["method"], raw_file={"path": "raw"}, checked_records=[])
    metadata = dict(search_recovery=True, participant="qfo_corrected_orthomcl_recovered",
                    query_coverage={"failed_queries": ["failed_fixture"]}, group_coverage={"fixture": True},
                    group_audit={"path": "group_audit"}, pair_semantics="cross_species_final_group_cliques")
    if recovered:
        audited.update(deepcopy(metadata))
    inputs = {"comparison": comparison, "baseline": baseline, "admission": report, "conversion": {},
              "replay": {"status": "corrected_checked_replay_admitted"}}
    # Admission and file integrity have dedicated tests; exercise assembly plus the real count validator.
    monkeypatch.setattr(module, "read_frozen", lambda path, digest: inputs[str(path)])
    monkeypatch.setattr(module, "record", lambda path: {"path": str(path)})
    monkeypatch.setattr(module, "check", lambda item: None)
    monkeypatch.setattr(module, "read_raw", lambda *args: ({}, [None] * counts["reference_relation_count"], {}))
    monkeypatch.setattr(module, "extract", lambda *args: {"key": target["method"]})
    monkeypatch.setattr(module, "audit_comparator", lambda *args: audited)
    result = module.audit(Path("comparison"), "digest", Path("baseline"))
    row = result["methods"][-1]
    assert row["families"] == target["families"]
    assert row["aggregate"] == target["aggregate"]
    assert result["publication_ready"] is False
    assert len(result["methods"]) == 8
    assert all(row["status"] == "not_admitted" for row in result["methods"][:-1])
    if recovered:
        assert {key: row[key] for key in metadata} == metadata
    else:
        assert not set(metadata).intersection(row)

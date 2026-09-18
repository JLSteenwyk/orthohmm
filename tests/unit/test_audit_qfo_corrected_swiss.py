from pathlib import Path

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

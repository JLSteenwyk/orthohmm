import copy
import json

import pytest

from benchmark_tools.audit_private_timing_advisories import audit, compare


def snapshot():
    alert = dict(number=1, state="open", ghsa_id="GHSA-test", vulnerable_version_range="<12.3.0",
                 dependency=dict(package=dict(name="Pillow", ecosystem="pip")))
    return dict(repository="JLSteenwyk/orthohmm", state_filter="open", alerts=[alert])


def test_duplicate_advisory_across_manifests_counted_once():
    old = snapshot()
    new = copy.deepcopy(old)
    new["alerts"].append(dict(new["alerts"][0], number=2))
    result = compare({"pillow": "12.2.0"}, [{"name": "Pillow", "version": "12.2.0"}], new, old)
    assert result["affected_manifest_alerts"] == 2
    assert result["affected_unique_advisories"] == 1
    assert result["added_alert_numbers"] == [2]


def test_fixed_and_absent_are_distinct():
    fixed = compare({"pillow": "12.3.0"}, [{"name": "Pillow", "version": "12.3.0"}], snapshot(), snapshot())
    absent = compare({}, [], snapshot(), snapshot())
    assert fixed["comparisons"][0]["status"] == "outside_reported_range"
    assert absent["comparisons"][0]["status"] == "not_in_lock"


def test_inventory_drift_rejected():
    with pytest.raises(ValueError, match="inventory"):
        compare({"pillow": "12.2.0"}, [], snapshot(), snapshot())


@pytest.mark.parametrize("status", ["private_packaging_historical_payload_aligned",
                                   "private_timing_environment_candidate_installed"])
def test_installed_receipt_audit(tmp_path, monkeypatch, status):
    python = tmp_path / "python"
    python.write_bytes(b"test interpreter identity")
    candidate = tmp_path / "candidate.json"
    candidate.write_text(json.dumps(dict(status=status, selected={"pillow": "12.3.0"})))
    alerts = tmp_path / "alerts.json"
    alerts.write_text(json.dumps(snapshot()))
    monkeypatch.setattr("benchmark_tools.audit_private_timing_advisories.subprocess.check_output",
                        lambda *a, **kw: '[{"name": "Pillow", "version": "12.3.0"}]')
    result = audit(python, candidate, alerts, alerts)
    assert result["affected_unique_advisories"] == 0
    assert result["scientific_execution_authorized"] is False
    assert result["comprehensive_security_clearance"] is False


@pytest.mark.parametrize("change", ["empty", "duplicate", "closed", "repository", "ecosystem"])
def test_invalid_snapshot_rejected(change):
    new = snapshot()
    if change == "empty":
        new["alerts"] = []
    elif change == "duplicate":
        new["alerts"] *= 2
    elif change == "closed":
        new["alerts"][0]["state"] = "dismissed"
    elif change == "repository":
        new["repository"] = "different/repo"
    else:
        new["alerts"][0]["dependency"]["package"]["ecosystem"] = "npm"
    with pytest.raises(ValueError, match="Require"):
        compare({}, [], new, snapshot())

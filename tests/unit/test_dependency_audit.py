from benchmark_tools.audit_dependency_lock import evaluate
from benchmark_tools.snapshot_dependency_alerts import summarize


def alert():
    return {"number": 1, "state": "open", "dependency": {"package": {"name": "Jinja2"}},
            "security_advisory": {"ghsa_id": "GHSA-test", "cve_id": None, "severity": "high", "summary": "Test",
                                  "description": "not retained", "references": [{"url": "https://example.org"}]},
            "security_vulnerability": {"vulnerable_version_range": "<=3.1.5", "first_patched_version": {"identifier": "3.1.6"}},
            "html_url": "https://example.org/1", "unexpected_secret": "never retain"}


def test_snapshot_uses_allowlist_and_omits_long_description():
    result = summarize(alert())
    assert "unexpected_secret" not in result and "description" not in result
    assert result["ghsa_id"] == "GHSA-test"


def test_every_lock_branch_checked_and_names_normalized():
    lock = {"package": [{"name": "jinja2", "version": "3.1.4"}, {"name": "jinja2", "version": "3.1.6"}]}
    row = evaluate(lock, [summarize(alert())])[0]
    assert row["affected_locked_versions"] == ["3.1.4"] and row["status"] == "affected"


def test_patched_version_is_outside_range():
    lock = {"package": [{"name": "jinja2", "version": "3.1.6"}]}
    assert evaluate(lock, [summarize(alert())])[0]["status"] == "outside_reported_range"


def test_absent_dependency_is_explicit_not_a_patched_claim():
    assert evaluate({"package": []}, [summarize(alert())])[0]["status"] == "not_in_lock"

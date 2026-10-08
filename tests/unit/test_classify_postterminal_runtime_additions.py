import copy
from datetime import datetime, timedelta, timezone

import pytest

from benchmark_tools.classify_postterminal_runtime_additions import classify


END = datetime(2026, 10, 8, 12, 56, tzinfo=timezone.utc)
CHECK = END + timedelta(seconds=30)
INSTALL = END + timedelta(hours=1)
ADDED = [dict(path="/usr/bin/lftp", mode=0o755, kind="file", bytes=1700488, sha256="a" * 64)]


def comparison():
    return dict(equal=False, expected_records=39208, observed_records=39209,
                added=copy.deepcopy(ADDED), missing=[], changed=[], metadata_differences={})


def test_exact_additions_classified_without_admission():
    result = classify(comparison(), ADDED, END, CHECK, INSTALL)
    assert result["status"] == "postterminal_additions_classified_not_admitted"
    assert result["original_inventory_equality"] is False
    for key in ("terminal_reviewed", "next_identity_authorized", "accuracy_evaluated",
                "publication_ready", "historical_failure_cause_established", "automatic_retry"):
        assert result[key] is False


@pytest.mark.parametrize("key,value", [
    ("equal", True), ("missing", [{"path": "/usr/bin/python"}]),
    ("changed", [{"path": "/usr/bin/python"}]), ("metadata_differences", {"roots": "changed"}),
    ("added", []), ("observed_records", 39210), ("expected_records", True),
])
def test_any_unapproved_difference_refused(key, value):
    case = comparison()
    case[key] = value
    with pytest.raises(ValueError):
        classify(case, ADDED, END, CHECK, INSTALL)


@pytest.mark.parametrize("end,check,installed", [
    (END.replace(tzinfo=None), CHECK, INSTALL), (END, CHECK.replace(tzinfo=None), INSTALL),
    (END, CHECK, INSTALL.replace(tzinfo=None)), (END, CHECK, END),
    (END, CHECK, CHECK), (CHECK, END, INSTALL),
])
def test_missing_or_non_postterminal_chronology_refused(end, check, installed):
    with pytest.raises(ValueError):
        classify(comparison(), ADDED, end, check, installed)


@pytest.mark.parametrize("key,value", [
    ("kind", "symlink"), ("path", "relative/lftp"), ("mode", True), ("bytes", 0),
    ("sha256", "A" * 64), ("sha256", "a" * 63),
])
def test_bad_exact_addition_descriptor_refused(key, value):
    added = copy.deepcopy(ADDED)
    added[0][key] = value
    case = comparison()
    case["added"] = added
    with pytest.raises(ValueError):
        classify(case, added, END, CHECK, INSTALL)


def test_duplicate_additional_paths_refused():
    added = ADDED * 2
    case = comparison()
    case.update(added=added, observed_records=39210)
    with pytest.raises(ValueError, match="sorted unique"):
        classify(case, added, END, CHECK, INSTALL)

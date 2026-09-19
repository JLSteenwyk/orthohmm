import gzip
import json
from pathlib import Path
from statistics import median

import pytest

from benchmark_tools import describe_root_context_controls as module

RESULTS = Path(module.__file__).parent / "results"
AUDIT = RESULTS / "root_context_audit_22019_20260919.json.gz"
SHA = "06bdeb7334454c4e616fbf783334d33bd5c72b3b3e73aafd7e452f87d24f639f"


@pytest.fixture
def audit():
    return json.loads(gzip.decompress(AUDIT.read_bytes()))


def test_signed_distributions_preserve_negatives():
    assert module.distributions([dict(x=-3), dict(x=-1), dict(x=2)]) == dict(
        x=dict(n=3, minimum=-3, median=-1, maximum=2))


@pytest.mark.parametrize("rows", [[], [dict(x=1), dict(y=1)], [dict(x=float("nan"))], [dict(x=True)]])
def test_invalid_distributions_rejected(rows):
    with pytest.raises(ValueError):
        module.distributions(rows)


def test_actual_stopped_panel_retains_all_missing_conditions(audit):
    result = module.describe(audit)
    assert len(result["trials"]) == 12
    assert [row["status"] for row in result["trials"][:4]] == ["validated"]*3 + ["retained_failure"]
    assert all(row["common"] is None for row in result["trials"][3:])
    assert result["blocks"][0]["common_interval_median_differences"]["user-contended_minus_steady"] is None
    assert all(value is None for block in result["blocks"][1:] for value in block["common_interval_median_differences"].values())
    assert result["scientific_timings_admitted"] is False


def test_actual_distributions_match_raw_comparisons(audit):
    result = module.describe(audit)
    for row, described in zip(audit["trials"][:3], result["trials"]):
        replay = row["replay"]
        comparisons = replay["measurement"]["context"]["intervals"]
        selected = replay["trial"]["common_intervals"]
        values = [comparisons[i]["root_minus_three_named_children_cpu_usec"] for i in selected]
        stats = described["common"]["distributions"]["root_minus_three_named_children_cpu_usec"]
        assert stats == dict(n=len(values), minimum=min(values), median=median(values), maximum=max(values))
        native = replay["measurement"]["lineage"]
        assert described["common"]["original_flags"] == len(set(selected) & set(native["original_flagged_intervals"]))
        assert described["all"]["narrow_flags"] == len(native["narrow_flagged_intervals"])


@pytest.mark.parametrize("fault", ["order", "count", "common", "intervals"])
def test_inconsistent_audit_rejected(audit, fault):
    if fault == "order":
        audit["trials"].reverse()
    elif fault == "count":
        audit["validated_trials"] = 12
    elif fault == "common":
        audit["trials"][0]["replay"]["trial"]["common_intervals"] = [99999]
    else:
        audit["trials"][0]["replay"]["measurement"]["context"]["intervals"].pop()
    with pytest.raises(ValueError):
        module.describe(audit)


def test_pinned_report_and_retained_description():
    result = module.report(AUDIT, SHA)
    retained = json.loads((RESULTS / "root_context_description_22019_20260919.json").read_text())
    for key in ("audit", "source"):
        result[key].pop("path")
        retained[key].pop("path")
    assert result == retained
    with pytest.raises(ValueError, match="pinned"):
        module.report(AUDIT, "0"*64)

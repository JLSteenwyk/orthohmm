import pytest

from benchmark_tools.probe_narrow_band_rescue import BANDS, summarize


def fixture():
    return {"pairs_per_band": 2, "frozen_runtime_modified": False, "bands": [
        {"band": band, "baseline": [1, 3], "scalar": [1, 2], "jit": [1, 2], "variant": [1, 2],
         "baseline_scalar_mismatches": 1, "variant_scalar_mismatches": 0,
         "variant_jit_mismatches": 0, "changed_indices": [1]} for band in BANDS]}


def test_summary_checks_raw_scores():
    result = summarize(fixture())
    assert result["variant_matches_scalar_and_jit"] and result["baseline_discrepancy_observed"]
    assert not result["frozen_runtime_modified"] and not result["benchmark_admitted"]


@pytest.mark.parametrize("change", ["missing_band", "short_scores", "float_scores", "false_count", "false_changes"])
def test_corrupt_diagnostic_rejected(change):
    report = fixture()
    row = report["bands"][0]
    if change == "missing_band":
        report["bands"].pop()
    elif change == "short_scores":
        row["jit"] = [1]
    elif change == "float_scores":
        row["variant"] = [1., 2.]
    elif change == "false_count":
        row["variant_scalar_mismatches"] = 1
    else:
        row["changed_indices"] = []
    with pytest.raises(ValueError):
        summarize(report)

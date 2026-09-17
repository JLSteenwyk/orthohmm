import pytest

from benchmark_tools.probe_narrow_band_rescue import BANDS, summarize, boundary_fixtures, orderings


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


def test_boundary_fixture_contains_cutoff_neighbors_and_empty_targets():
    queries, targets, pairs = boundary_fixtures()
    assert {49, 50, 51} <= {len(q) for q in queries}
    assert {0, 49, 50, 51} <= {len(t) for t in targets}
    assert pairs.shape == (495, 2)
    assert len(set(map(tuple, pairs))) == 495
    assert all(sorted(order.tolist()) == list(range(495)) for order in orderings(495).values())


@pytest.mark.parametrize("change", [None, "missing", "false_score", "false_indices"])
def test_order_summary_recomputes_scores(change):
    report = fixture()
    report["fixture_kind"] = "boundary_lengths"
    report["order_checks"] = [
        {"band": band, "order": name, "threads": threads, "mismatches": 0,
         "mismatched_original_indices": [], "scores": [[1, 2][i] for i in order]}
        for band in BANDS for name, order in orderings(2).items() for threads in (1, 4)]
    if change == "missing":
        report["order_checks"].pop()
    elif change == "false_score":
        report["order_checks"][0]["scores"][0] = 99
    elif change == "false_indices":
        report["order_checks"][0]["mismatched_original_indices"] = [0]
    if change:
        with pytest.raises(ValueError):
            summarize(report)
    else:
        assert len(summarize(report)["order_checks"]) == 30

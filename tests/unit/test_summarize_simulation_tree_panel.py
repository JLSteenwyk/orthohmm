from copy import deepcopy

import pytest

from benchmark_tools import summarize_simulation_tree_panel as summary


def score(tp=5, fp=5, fn=5, scale=1):
    tp, fp, fn = tp * scale, fp * scale, fn * scale
    return {"tp": tp, "fp": fp, "fn": fn, "input_genes": 10000, "eligible_true_pairs": tp + fn,
            "f1": 2 * tp / (2 * tp + fp + fn), "precision": tp / (tp + fp), "recall": tp / (tp + fn),
            "undefined_ratios": []}


def rows():
    return [{"condition": c, "seed": s, "method": m, "arm": a, "status": "complete", "truth_sha256": "a" * 64,
             "score": score()} for c in summary.CONDITIONS for s in summary.SEEDS for m in summary.METHODS for a in summary.ARMS]


@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "unknown", "pending", "truth", "counts", "imputed"])
def test_full_panel_contract(problem):
    records = rows()
    if problem == "missing":
        records.pop()
    elif problem == "duplicate":
        records[-1] = deepcopy(records[0])
    elif problem == "unknown":
        records[0]["arm"] = "selected_best"
    elif problem == "pending":
        records[0]["status"] = "running"
    elif problem == "truth":
        records[0]["truth_sha256"] = "b" * 64
    elif problem == "counts":
        records[0]["score"]["tp"] += 1
    elif problem == "imputed":
        records[0].update(status="failed", reason="fixture")
    if problem:
        with pytest.raises(ValueError):
            summary.index_records(records)
    else:
        assert len(summary.index_records(records)) == 560


def test_arithmetic_seed_mean_not_pooled_pair_counts():
    result = summary.paired_summary({1: score(9, 1, 1), 2: score(2, 8, 8, 100)},
                                    {1: score(), 2: score(scale=100)})
    assert result["metrics"]["f1"]["difference_percentage_points"] == pytest.approx(5.)
    assert result["metrics"]["f1"]["target_mean"] == pytest.approx(.55)
    nominal = result["metrics"]["f1"]["paired_95_percent_ci"]
    adjusted = result["metrics"]["f1"]["bonferroni_126_ci"]
    assert adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1]


def test_fixed_seed_and_sign_symmetry():
    a, b = {i: score(i, 10-i, 10-i) for i in range(1, 10)}, {i: score() for i in range(1, 10)}
    first = summary.paired_summary(a, b)
    assert first == summary.paired_summary(a, b)
    reverse = summary.paired_summary(b, a)
    for metric in summary.METRICS:
        x, y = first["metrics"][metric], reverse["metrics"][metric]
        assert x["difference_percentage_points"] == pytest.approx(-y["difference_percentage_points"])
        assert x["bonferroni_126_ci"] == pytest.approx([-v for v in reversed(y["bonferroni_126_ci"])])


def test_no_and_single_pairs_do_not_fabricate_intervals():
    empty = summary.paired_summary({}, {})
    assert empty["status"] == "no_complete_pairs" and len(empty["metrics"]) == 3
    one = summary.paired_summary({1: score()}, {1: score()})
    assert one["status"] == "insufficient_seeds"
    assert all(m["bonferroni_126_ci"] is None for m in one["metrics"].values())


def test_mismatched_pairs_rejected():
    with pytest.raises(ValueError):
        summary.paired_summary({1: score()}, {2: score()})


def test_all_126_slots_retained_even_when_every_run_failed():
    records = rows()
    for row in records:
        row.update(status="failed", reason="fixture")
        del row["score"]
    result = summary.summarize(records)
    assert len(result["contrasts"]) == 42 and result["bootstrap"]["multiplicity"] == 126
    assert sum(len(c["metrics"]) for c in result["contrasts"]) == 126
    assert all(c["paired_seed_count"] == 0 and len(c["excluded_seeds"]) == 10 for c in result["contrasts"])


def test_failure_excludes_only_affected_pair_and_never_changes_multiplicity():
    records = rows()
    records[0].update(status="failed", reason="fixture")
    del records[0]["score"]
    result = summary.summarize(records)
    first = result["contrasts"][0]
    assert first["paired_seed_count"] == 9 and len(first["excluded_seeds"]) == 1
    assert result["contrasts"][1]["paired_seed_count"] == 10
    assert result["bootstrap"]["multiplicity"] == 126

import copy

import numpy as np
import pytest

from benchmark_tools.summarize_simulation_panel import CONDITIONS, METHODS, SEEDS, paired_seed_summary, summarize


def score(tp=1, fp=0, fn=0):
    return {"tp": tp, "fp": fp, "fn": fn, "input_genes": 20, "eligible_true_pairs": tp + fn,
            "f1": 2 * tp / (2 * tp + fp + fn), "precision": tp / (tp + fp),
            "recall": tp / (tp + fn), "undefined_ratios": []}


def records():
    return [{"condition": c, "seed": s, "method": m, "status": "complete",
             "truth_sha256": "a" * 64, "score": score()} for c in CONDITIONS for s in SEEDS for m in METHODS]


def test_direct_paired_seed_resampling_and_not_pooled_counts():
    a, b = {1: score(1, 1), 2: score(100, 0)}, {1: score(1), 2: score(100)}
    result = paired_seed_summary(a, b)
    assert result["metrics"]["f1"]["comparator_mean"] == pytest.approx((2 / 3 + 1) / 2)
    assert result["metrics"]["f1"]["comparator_mean"] != pytest.approx(202 / 203)
    rng = np.random.Generator(np.random.PCG64(20261031))
    weights = rng.multinomial(2, [0.5, 0.5], size=20000)
    draws = weights[:, 0] * 100 * (2 / 3 - 1) / 2
    assert result["metrics"]["f1"]["paired_95_percent_ci"] == pytest.approx(np.quantile(draws, [.025, .975]))
    assert result["metrics"]["f1"]["bonferroni_14_ci"] == pytest.approx(np.quantile(draws, [.025 / 14, 1 - .025 / 14]))
    assert "bonferroni_14_ci" not in result["metrics"]["precision"]


def test_full_panel_identity_and_diagnostic_exclusion():
    report = summarize(records())
    for block in report["conditions"].values():
        assert set(block["contrasts"]) == set(METHODS[:2])
        for contrast in block["contrasts"].values():
            assert len(contrast["included_seeds"]) == 10
            assert contrast["metrics"]["f1"]["bonferroni_14_ci"] == [0, 0]


def test_failures_are_explicit_conditional_exclusions():
    rows = records()
    first = rows[0]
    first.update(status="failed", reason="test process failure")
    del first["score"]
    report = summarize(rows)
    contrast = report["conditions"]["baseline"]["contrasts"][METHODS[0]]
    assert contrast["included_seeds"] == list(SEEDS[1:])
    assert contrast["excluded_seeds"][0]["seed"] == SEEDS[0]
    assert contrast["conditional_on_success"] is True
    assert report["conditions"]["baseline"]["methods"][METHODS[0]]["failure_fraction_of_planned"] == .1


@pytest.mark.parametrize("fault", ["missing", "duplicate", "wrong_truth", "bad_score", "failure_with_score", "pending", "missing_reason"])
def test_invalid_panel_rejected(fault):
    rows = records()
    if fault == "missing":
        rows.pop()
    elif fault == "duplicate":
        rows.append(copy.deepcopy(rows[0]))
    elif fault == "wrong_truth":
        rows[0]["truth_sha256"] = "b" * 64
    elif fault == "bad_score":
        rows[0]["score"]["f1"] = .5
    elif fault == "failure_with_score":
        rows[0].update(status="failed", reason="test")
    elif fault == "pending":
        rows[0]["status"] = "pending"
    else:
        rows[0]["status"] = "failed"
        del rows[0]["score"]
    with pytest.raises(ValueError):
        summarize(rows)


def test_too_few_seeds_do_not_produce_spurious_intervals():
    assert paired_seed_summary({}, {})["status"] == "no_complete_pairs"
    one = paired_seed_summary({1: score()}, {1: score()})
    assert one["status"] == "insufficient_seeds"
    assert "paired_95_percent_ci" not in one["metrics"]["f1"]

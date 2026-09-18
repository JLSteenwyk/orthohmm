import copy

import numpy as np
import pytest

from benchmark_tools.bootstrap_corrected_swiss_strata import (
    CONTRASTS, METHODS, METRICS, bootstrap, validated_values,
)


def fixture():
    membership = {f"f{i}": "lower" if i < 9 else "higher" for i in range(18)}
    counts = {}
    for k, method in enumerate(METHODS):
        counts[method] = {}
        for i, family in enumerate(membership):
            tp, fp = (i + k) % 9, (2 * i + 3 * k) % 11
            counts[method][family] = dict(
                counts_without_prior=dict(TP=tp, FN=10-tp, FP=fp, TN=12-fp),
                represented_genes=[f"{family}g{j}" for j in range(8)])
    return counts, membership


def direct_scores(rows):
    # Independently reproduce native half-count-plus-one smoothing and macro F1.
    p = np.mean([(r["TP"] / 2 + 1) / ((r["TP"] + r["FP"]) / 2 + 2) for r in rows])
    r = np.mean([(r["TP"] / 2 + 1) / ((r["TP"] + r["FN"]) / 2 + 2) for r in rows])
    return np.array([2*p*r/(p+r), p, r])


def test_all_27_endpoints_match_independent_repeated_family_enumeration():
    counts, membership = fixture()
    report = bootstrap(counts, membership, replicates=200, seed=42)
    rng = np.random.Generator(np.random.PCG64(42))
    point, draws = {}, {}
    for name in ("lower", "higher"):
        families = [f for f in membership if membership[f] == name]
        weights = rng.multinomial(9, np.full(9, 1/9), size=200)
        point[name], draws[name] = [], []
        for method in METHODS:
            rows = [counts[method][f]["counts_without_prior"] for f in families]
            point[name].append(direct_scores(rows))
            draws[name].append([direct_scores([r for r, n in zip(rows, w) for _ in range(n)]) for w in weights])
        point[name], draws[name] = np.array(point[name]), np.array(draws[name])
        for i, method in enumerate(METHODS):
            assert list(report["bins"][name]["point_estimates"][method].values()) == pytest.approx(point[name][i])
    checked = 0
    for name in ("lower", "higher", "interaction"):
        comparisons = report["interactions"] if name == "interaction" else report["bins"][name]["comparisons"]
        for row, (candidate, reference) in zip(comparisons, CONTRASTS):
            if name == "interaction":
                expected = (draws["higher"][candidate] - draws["higher"][reference]
                            - (draws["lower"][candidate] - draws["lower"][reference]))
                observed = (point["higher"][candidate] - point["higher"][reference]
                            - (point["lower"][candidate] - point["lower"][reference]))
            else:
                expected = draws[name][candidate] - draws[name][reference]
                observed = point[name][candidate] - point[name][reference]
            for j, metric in enumerate(METRICS):
                result = row["metrics"][metric]
                assert result["difference"] == pytest.approx(observed[j])
                assert result["paired_percentile_ci"] == pytest.approx(np.quantile(expected[:, j], [.025, .975]))
                assert result["bonferroni_percentile_ci"] == pytest.approx(np.quantile(expected[:, j], [.05/54, 1-.05/54]))
                if name != "interaction":
                    assert sum(result[k] for k in ("family_wins", "family_ties", "family_losses")) == 9
                checked += 1
    assert checked == report["multiplicity_endpoints"] == 27
    assert report["scientific_inputs_admitted"] is False


@pytest.mark.parametrize("lower,missing", [(4, 0), (0, 0), (5, 8), (0, 18)])
def test_small_empty_and_missing_bins_never_get_intervals(lower, missing):
    counts, membership = fixture()
    for i, f in enumerate(membership):
        membership[f] = "lower" if i < lower else "missing" if i >= 18-missing else "higher"
    report = bootstrap(counts, membership, replicates=100)
    assert report["multiplicity_endpoints"] == 27
    for name, row in report["bins"].items():
        eligible = name != "missing" and len(row["families"]) >= 5
        assert row["interval_eligible"] == eligible
        for comparison in row["comparisons"]:
            for result in comparison["metrics"].values():
                assert (result["paired_percentile_ci"] is not None) == eligible
                assert (result["difference"] is not None) == bool(row["families"])
    for comparison in report["interactions"]:
        for result in comparison["metrics"].values():
            eligible = all(report["bins"][name]["interval_eligible"] for name in ("lower", "higher"))
            assert (result["paired_percentile_ci"] is not None) == eligible


def test_identical_methods_have_zero_contrasts_and_interactions():
    counts, membership = fixture()
    for method in METHODS[1:]:
        counts[method] = copy.deepcopy(counts[METHODS[0]])
    report = bootstrap(counts, membership, replicates=100)
    for name in ("lower", "higher"):
        for comparison in report["bins"][name]["comparisons"]:
            for result in comparison["metrics"].values():
                assert result["difference"] == 0
                assert result["bonferroni_percentile_ci"] == [0, 0]
                assert result["family_ties"] == 9
    for comparison in report["interactions"]:
        for result in comparison["metrics"].values():
            assert result["difference"] == 0
            assert result["bonferroni_percentile_ci"] == [0, 0]


@pytest.mark.parametrize("change", ["method", "family", "bin", "negative", "bool", "truth", "members", "overlap"])
def test_invalid_counts_or_membership_rejected(change):
    counts, membership = fixture()
    row = counts[METHODS[1]]["f0"]
    if change == "method":
        counts.pop(METHODS[0])
    elif change == "family":
        membership.pop("f0")
    elif change == "bin":
        membership["f0"] = "unknown"
    elif change in ("negative", "bool"):
        row["counts_without_prior"]["TP"] = -1 if change == "negative" else True
    elif change == "truth":
        row["counts_without_prior"]["TN"] += 1
    elif change == "members":
        row["represented_genes"][0] = "foreign"
    else:
        row["represented_genes"][0] = "f1g0"
    with pytest.raises(ValueError):
        validated_values(counts, membership)


@pytest.mark.parametrize("kwargs", [dict(replicates=True), dict(replicates=99), dict(seed=-1), dict(seed=True)])
def test_invalid_resampling_controls_rejected(kwargs):
    counts, membership = fixture()
    with pytest.raises(ValueError):
        bootstrap(counts, membership, **kwargs)

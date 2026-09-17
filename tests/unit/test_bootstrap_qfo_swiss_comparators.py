import copy

import numpy as np
import pytest

from benchmark_tools.audit_qfo_swiss_comparators import METHODS
from benchmark_tools.audit_qfo_swiss_counts import statistics
from benchmark_tools.bootstrap_qfo_swiss_comparators import CONTRASTS, METRICS, bootstrap, validated_values


def fixture():
    families = [f"f{i}" for i in range(18)]
    methods = []
    for k, name in enumerate(METHODS):
        rows = []
        for i, family in enumerate(families):
            tp, fp = 1 + (i + k) % 5, 1 + (i * 2 + k) % 5
            counts = dict(TP=tp, FN=6 - tp, FP=fp, TN=6 - fp)
            rows.append(dict(family=family, counts_without_prior=counts,
                             statistics_with_prior=statistics(counts),
                             represented_genes=[f"{family}g{g}" for g in range(6)]))
        p = sum(r["statistics_with_prior"]["PPV"] for r in rows) / 18
        r = sum(r["statistics_with_prior"]["TPR"] for r in rows) / 18
        methods.append(dict(method=name, families=rows, aggregate=dict(PPV=p, TPR=r, F1=2*p*r/(p+r))))
    return dict(status="eight_comparator_swiss_family_counts_verified", shared_represented_genes={},
                families=families, methods=methods, reference_relation_count=216)


def test_direct_repeated_family_enumeration_matches_all_intervals():
    report = fixture()
    result = bootstrap(report, replicates=200, seed=42)
    weights = np.random.Generator(np.random.PCG64(42)).multinomial(18, np.full(18, 1/18), size=200)
    draws = []
    for method in report["methods"]:
        method_draws = []
        for draw in weights:
            selected = [row for row, n in zip(method["families"], draw) for _ in range(n)]
            p = sum(row["statistics_with_prior"]["PPV"] for row in selected) / 18
            r = sum(row["statistics_with_prior"]["TPR"] for row in selected) / 18
            method_draws.append([2*p*r/(p+r), p, r])
        draws.append(method_draws)
    draws = np.asarray(draws)
    assert result["multiplicity_endpoints"] == 24
    for row, (candidate, reference) in zip(result["comparisons"], CONTRASTS):
        differences = draws[candidate] - draws[reference]
        for j, metric in enumerate(METRICS):
            m = row["metrics"][metric]
            assert m["paired_percentile_ci"] == pytest.approx(np.quantile(differences[:, j], [.025, .975]))
            assert m["bonferroni_percentile_ci"] == pytest.approx(np.quantile(differences[:, j], [.05/48, 1-.05/48]))
            assert sum(m[k] for k in ("family_wins", "family_ties", "family_losses")) == 18


def test_identical_methods_have_exact_zero_differences():
    report = fixture()
    for method in report["methods"]:
        method["families"] = copy.deepcopy(report["methods"][0]["families"])
        method["aggregate"] = report["methods"][0]["aggregate"].copy()
    result = bootstrap(report, replicates=100)
    for row in result["comparisons"]:
        for m in row["metrics"].values():
            assert m["difference"] == 0
            assert m["paired_percentile_ci"] == [0, 0]
            assert m["bonferroni_percentile_ci"] == [0, 0]
            assert m["family_ties"] == 18


@pytest.mark.parametrize("change", ["order", "negative", "bool", "overlap", "stat", "aggregate", "truth"])
def test_invalid_sufficient_statistics_rejected(change):
    report = fixture()
    row = report["methods"][1]["families"][0]
    if change == "order":
        report["methods"].reverse()
    elif change == "negative":
        row["counts_without_prior"]["TP"] = -1
    elif change == "bool":
        row["counts_without_prior"]["TP"] = True
    elif change == "overlap":
        row["represented_genes"][0] = "f1g0"
    elif change == "stat":
        row["statistics_with_prior"]["F1"] += .1
    elif change == "aggregate":
        report["methods"][1]["aggregate"]["PPV"] += .1
    else:
        row["counts_without_prior"]["TN"] += 1
    with pytest.raises(ValueError):
        validated_values(report)

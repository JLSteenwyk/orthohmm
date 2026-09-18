from copy import deepcopy

import numpy as np
import pytest

from benchmark_tools.audit_qfo_swiss_counts import REFERENCE_SHA, statistics
from benchmark_tools.bootstrap_qfo_factorial import CELLS, METRICS, contrasts, bootstrap, validated_values, markdown


def fixture():
    families = [f"family{i}" for i in range(18)]
    cells = []
    for index, label in enumerate(CELLS):
        rows = []
        for i, family in enumerate(families):
            tp, fp = 6 + i % 3 + index % 4, 4 + i % 4 + index // 2
            counts = {"TP": tp, "FN": 20 - tp, "FP": fp, "TN": 20 - fp}
            rows.append({"family": family, "represented_genes": [f"{family}_gene{j}" for j in range(12)],
                         "counts_without_prior": counts, "statistics_with_prior": statistics(counts)})
        p = sum(r["statistics_with_prior"]["PPV"] for r in rows) / 18
        r = sum(r["statistics_with_prior"]["TPR"] for r in rows) / 18
        cells.append({"cell": label, "families": rows, "aggregate": {"F1": 2 * p * r / (p + r), "PPV": p, "TPR": r}})
    return {"status": "qfo_factorial_swiss_counts_verified", "reference": {"sha256": REFERENCE_SHA},
            "shared_represented_genes": {}, "reference_relation_count": 720, "families": families, "cells": cells}


def test_contrast_inventory_and_orientation():
    rows = contrasts()
    assert len(rows) == 14 and len({r["name"] for r in rows}) == 14
    edges = set()
    for row in rows[:12]:
        weights = row["weights"]
        baseline, candidate = weights.index(-1), weights.index(1)
        assert sum(abs(w) for w in weights) == 2
        assert candidate - baseline == {"P": 4, "C": 2, "R": 1}[row["name"][0]]
        assert CELLS[candidate] == row["candidate"] and CELLS[baseline] == row["reference"]
        edges.add((baseline, candidate))
    assert len(edges) == 12
    assert rows[12]["weights"] == [1, -1, -1, 1, 0, 0, 0, 0]
    assert rows[13]["weights"] == [0, 0, 0, 0, 1, -1, -1, 1]


def direct_metrics(rows):
    p, r = [], []
    for row in rows:
        c = row["counts_without_prior"]
        tp, fp, fn = c["TP"] / 2 + 1, c["FP"] / 2 + 1, c["FN"] / 2 + 1
        p.append(tp / (tp + fp))
        r.append(tp / (tp + fn))
    p, r = sum(p) / len(p), sum(r) / len(r)
    return np.array([2 * p * r / (p + r), p, r])


def test_all_42_intervals_match_explicit_repeated_family_reaggregation():
    data = fixture()
    result = bootstrap(data, replicates=100, seed=57)
    weights = np.random.Generator(np.random.PCG64(57)).multinomial(18, np.full(18, 1 / 18), size=100)
    raw_draws = []
    for cell in data["cells"]:
        raw_draws.append([direct_metrics([row for row, count in zip(cell["families"], draw) for _ in range(count)]) for draw in weights])
    raw_draws = np.asarray(raw_draws)
    points = np.asarray([direct_metrics(cell["families"]) for cell in data["cells"]])
    for contrast in result["comparisons"]:
        # Explicit signed sums are independent of the implementation's tensor contraction.
        delta = sum(w * raw_draws[i] for i, w in enumerate(contrast["weights"]))
        point = sum(w * points[i] for i, w in enumerate(contrast["weights"]))
        for j, metric in enumerate(METRICS):
            actual = contrast["metrics"][metric]
            assert actual["difference"] == pytest.approx(point[j], abs=1e-14)
            for key, q in (("paired_percentile_ci", [.025, .975]), ("bonferroni_percentile_ci", [.05 / 84, 1 - .05 / 84])):
                assert actual[key] == pytest.approx(np.quantile(delta[:, j], q, method="linear"), abs=1e-14)
            assert sum(actual[k] for k in ("family_wins", "family_ties", "family_losses")) == 18
    assert result["multiplicity_endpoints"] == 42
    assert len(markdown(result).splitlines()) >= 42


def test_identical_cells_have_exact_zero_contrasts():
    data = fixture()
    for cell in data["cells"]:
        cell["families"] = deepcopy(data["cells"][0]["families"])
        cell["aggregate"] = deepcopy(data["cells"][0]["aggregate"])
    result = bootstrap(data, 100, 58)
    for contrast in result["comparisons"]:
        for metric in contrast["metrics"].values():
            assert metric["difference"] == 0
            assert metric["bonferroni_percentile_ci"] == [0, 0]
            assert metric["family_ties"] == 18


@pytest.mark.parametrize("mutation", ["missing_cell", "order", "missing_family", "reference", "overlap",
    "negative", "float", "truth", "members", "family_stat", "aggregate", "relations"])
def test_invalid_count_inputs_rejected(mutation):
    data = fixture()
    row = data["cells"][1]["families"][0]
    if mutation == "missing_cell":
        data["cells"].pop()
    elif mutation == "order":
        data["cells"].reverse()
    elif mutation == "missing_family":
        data["families"].pop()
    elif mutation == "reference":
        data["reference"]["sha256"] = "other"
    elif mutation == "overlap":
        data["cells"][1]["families"][1]["represented_genes"] = row["represented_genes"]
    elif mutation in ("negative", "float", "truth"):
        row["counts_without_prior"]["TP"] = {"negative": -1, "float": 7.0, "truth": 8}[mutation]
        row["statistics_with_prior"] = statistics(row["counts_without_prior"])
    elif mutation == "members":
        row["represented_genes"][0] = "changed_gene"
    elif mutation == "family_stat":
        row["statistics_with_prior"]["PPV"] += .01
    elif mutation == "aggregate":
        data["cells"][1]["aggregate"]["F1"] += .01
    else:
        data["reference_relation_count"] += 1
    with pytest.raises(ValueError):
        validated_values(data)

from copy import deepcopy

import numpy as np
import pytest

from benchmark_tools import bootstrap_qfo_sequence as module
from benchmark_tools.audit_qfo_swiss_counts import statistics
from tests.unit.test_bootstrap_qfo_factorial import fixture as factorial, direct_metrics


def fixture():
    data = factorial()
    variants = []
    for name, cell in zip(module.VARIANTS, data.pop("cells")[:3]):
        cell.pop("cell")
        cell["variant"] = name
        for index, row in enumerate(cell["families"]):
            row["counts_without_prior"]["TN"] += 558 if index < 17 else 559
            row["statistics_with_prior"] = statistics(row["counts_without_prior"])
        variants.append(cell)
    data.update(status="corrected_qfo_sequence_swiss_counts_verified", variants=variants,
                reference_relation_count=10765, publication_ready=False, uncertainty_admitted=False)
    return data


def test_explicit_resampling_and_all_six_endpoints():
    data = fixture()
    saved = deepcopy(data)
    result = module.bootstrap(data, replicates=100, seed=57)
    assert data == saved
    draws = np.random.Generator(np.random.PCG64(57)).multinomial(18, np.full(18, 1 / 18), size=100)
    scores = np.asarray([[direct_metrics([row for row, n in zip(v["families"], draw) for _ in range(n)])
                         for draw in draws] for v in data["variants"]])
    points = np.asarray([direct_metrics(v["families"]) for v in data["variants"]])
    for index, contrast in enumerate(result["comparisons"], 1):
        assert contrast["candidate"] == module.VARIANTS[index]
        assert contrast["reference"] == "p0_c0_r0"
        for j, metric in enumerate(module.METRICS):
            actual = contrast["metrics"][metric]
            assert actual["difference"] == pytest.approx(points[index, j] - points[0, j], abs=1e-14)
            for key, q in (("paired_percentile_ci", [.025, .975]),
                           ("bonferroni_percentile_ci", [.05 / 12, 1 - .05 / 12])):
                assert actual[key] == pytest.approx(np.quantile(scores[index, :, j] - scores[0, :, j], q), abs=1e-14)
            assert sum(actual[k] for k in ("family_wins", "family_ties", "family_losses")) == 18
    assert result["multiplicity_endpoints"] == 6
    assert result["uncertainty_admitted"] is False


def test_identical_and_zero_prediction_counts():
    data = fixture()
    rows = data["variants"][0]["families"]
    for row in rows:
        c = row["counts_without_prior"]
        c["FN"] += c["TP"]
        c["TN"] += c["FP"]
        c["TP"] = c["FP"] = 0
        row["statistics_with_prior"] = statistics(c)
    for variant in data["variants"]:
        variant["families"] = deepcopy(rows)
        variant["aggregate"] = dict(zip(module.METRICS, direct_metrics(rows)))
    result = module.bootstrap(data, 100, 58)
    for contrast in result["comparisons"]:
        for value in contrast["metrics"].values():
            assert value["difference"] == 0
            assert value["bonferroni_percentile_ci"] == [0, 0]
            assert value["family_ties"] == 18


@pytest.mark.parametrize("mutation", ["status", "admitted", "publication", "reference", "relations",
    "missing", "order", "family_order", "overlap", "negative", "float", "boolean", "truth", "members", "stat", "aggregate"])
def test_rejects_invalid_counts(mutation):
    data = fixture()
    row = data["variants"][1]["families"][0]
    if mutation == "status":
        data["status"] = "historical"
    elif mutation == "admitted":
        data["uncertainty_admitted"] = True
    elif mutation == "publication":
        data["publication_ready"] = True
    elif mutation == "reference":
        data["reference"]["sha256"] = "other"
    elif mutation == "relations":
        data["reference_relation_count"] -= 1
    elif mutation == "missing":
        data["variants"].pop()
    elif mutation == "order":
        data["variants"].reverse()
    elif mutation == "family_order":
        data["variants"][1]["families"].reverse()
    elif mutation == "overlap":
        data["variants"][1]["families"][1]["represented_genes"] = row["represented_genes"]
    elif mutation in ("negative", "float", "boolean", "truth"):
        row["counts_without_prior"]["TP"] = {"negative": -1, "float": 1.5, "boolean": True, "truth": 11}[mutation]
        row["statistics_with_prior"] = statistics(row["counts_without_prior"])
    elif mutation == "members":
        row["represented_genes"][0] = "changed"
    elif mutation == "stat":
        row["statistics_with_prior"]["PPV"] = np.nan
    else:
        data["variants"][1]["aggregate"]["F1"] += .01
    with pytest.raises(ValueError):
        module.bootstrap(data, 100, 57)


@pytest.mark.parametrize("replicates,seed", [(True, 1), (99, 1), (100, -1), (100, False)])
def test_invalid_controls(replicates, seed):
    with pytest.raises(ValueError):
        module.bootstrap(fixture(), replicates, seed)

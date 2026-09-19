from copy import deepcopy

import numpy as np
import pytest

from benchmark_tools import bootstrap_qfo_parameter_neighborhood as module


def fixture():
    families = [f"family{i}" for i in range(18)]
    arms = []
    for index, label in enumerate(module.ARMS):
        rows = []
        for i, family in enumerate(families):
            tp, fp = 6 + i % 3 + index, 4 + i % 4 + index // 2
            counts = {"TP": tp, "FN": 20 - tp, "FP": fp, "TN": 20 - fp}
            rows.append({"family": family, "represented_genes": [f"{family}_g{j}" for j in range(12)],
                "counts_without_prior": counts, "statistics_with_prior": module.statistics(counts)})
        p = sum(r["statistics_with_prior"]["PPV"] for r in rows) / 18
        r = sum(r["statistics_with_prior"]["TPR"] for r in rows) / 18
        arms.append({"arm": label, "status": "counts_verified", "families": rows,
                     "aggregate": {"F1": 2 * p * r / (p + r), "PPV": p, "TPR": r}})
    return {"status": "corrected_qfo_parameter_swiss_counts_verified", "publication_ready": False,
        "uncertainty_admitted": False, "reference": {"sha256": module.REFERENCE_SHA},
        "shared_represented_genes": {}, "reference_relation_count": 720, "families": families, "arms": arms}


def direct(rows):
    precision, recall = [], []
    for row in rows:
        c = row["counts_without_prior"]
        tp, fp, fn = c["TP"] / 2 + 1, c["FP"] / 2 + 1, c["FN"] / 2 + 1
        precision.append(tp / (tp + fp))
        recall.append(tp / (tp + fn))
    p, r = sum(precision) / len(rows), sum(recall) / len(rows)
    return np.array([2 * p * r / (p + r), p, r])


def test_all_18_intervals_match_explicit_family_reaggregation():
    data = fixture()
    result = module.calculate(data, replicates=100, seed=57)
    weights = np.random.Generator(np.random.PCG64(57)).multinomial(18, np.full(18, 1 / 18), size=100)
    draws = np.asarray([[direct([row for row, count in zip(arm["families"], draw) for _ in range(count)])
                         for draw in weights] for arm in data["arms"]])
    points = np.asarray([direct(arm["families"]) for arm in data["arms"]])
    assert result["multiplicity_endpoints"] == 18
    assert result["protocol_controls_match"] is False
    for index, contrast in enumerate(result["comparisons"], 1):
        assert contrast["candidate"] == module.ARMS[index]
        for j, metric in enumerate(module.METRICS):
            actual = contrast["metrics"][metric]
            delta = draws[index, :, j] - draws[0, :, j]
            assert actual["difference"] == pytest.approx(points[index, j] - points[0, j], abs=1e-14)
            for key, q in (("paired_percentile_ci", [.025, .975]),
                           ("bonferroni_percentile_ci", [.05 / 36, 1 - .05 / 36])):
                assert actual[key] == pytest.approx(np.quantile(delta, q, method="linear"), abs=1e-14)
            assert sum(actual[k] for k in ("family_wins", "family_ties", "family_losses")) == 18


def test_missing_variants_do_not_change_shared_draws_or_multiplicity():
    data = fixture()
    complete = module.calculate(data, replicates=100, seed=58)
    for i in (1, 2):
        data["arms"][i] = {"arm": module.ARMS[i], "status": "not_admitted", "reason": "pending execution"}
    partial = module.calculate(data, replicates=100, seed=58)
    assert partial["multiplicity_endpoints"] == 18
    assert partial["comparisons"][2:] == complete["comparisons"][2:]
    for row in partial["comparisons"][:2]:
        assert row["status"] == "not_estimable" and row["metrics"] is None
    assert partial["point_estimates"]["cpm_low"] is None


def test_unavailable_baseline_prevents_all_paired_intervals():
    data = fixture()
    data["arms"][0] = {"arm": "control", "status": "not_admitted", "reason": "failed audit"}
    result = module.calculate(data, replicates=100)
    assert all(r["metrics"] is None and r["reason"] == "baseline_not_admitted" for r in result["comparisons"])
    assert result["point_estimates"]["norm_low"] is not None


def test_identical_arms_have_exact_zero_and_18_ties():
    data = fixture()
    for arm in data["arms"][1:]:
        arm["families"] = deepcopy(data["arms"][0]["families"])
        arm["aggregate"] = deepcopy(data["arms"][0]["aggregate"])
    result = module.calculate(data, replicates=100)
    for comparison in result["comparisons"]:
        for metric in comparison["metrics"].values():
            assert metric["difference"] == 0 and metric["bonferroni_percentile_ci"] == [0, 0]
            assert metric["family_ties"] == 18


@pytest.mark.parametrize("problem", ["status", "ready", "uncertainty", "reference", "order", "missing_arm",
    "missing_family", "overlap", "float", "negative", "truth", "members", "family_stat", "aggregate",
    "relations", "missing_reason", "imputed_missing"])
def test_invalid_counts_rejected(problem):
    data = fixture()
    row = data["arms"][1]["families"][0]
    if problem == "status":
        data["status"] = "unverified"
    elif problem == "ready":
        data["publication_ready"] = True
    elif problem == "uncertainty":
        data["uncertainty_admitted"] = True
    elif problem == "reference":
        data["reference"]["sha256"] = "other"
    elif problem == "order":
        data["arms"].reverse()
    elif problem == "missing_arm":
        data["arms"].pop()
    elif problem == "missing_family":
        data["families"].pop()
    elif problem == "overlap":
        data["arms"][1]["families"][1]["represented_genes"] = row["represented_genes"]
    elif problem in ("float", "negative", "truth"):
        row["counts_without_prior"]["TP"] = {"float": 7.0, "negative": -1, "truth": 8}[problem]
        row["statistics_with_prior"] = module.statistics(row["counts_without_prior"])
    elif problem == "members":
        row["represented_genes"][0] = "changed"
    elif problem == "family_stat":
        row["statistics_with_prior"]["PPV"] += .01
    elif problem == "aggregate":
        data["arms"][1]["aggregate"]["F1"] += .01
    elif problem == "relations":
        data["reference_relation_count"] += 1
    else:
        data["arms"][1] = {"arm": "cpm_low", "status": "not_admitted"}
        if problem == "imputed_missing":
            data["arms"][1].update(reason="failed", aggregate={"F1": 0})
    with pytest.raises(ValueError):
        module.calculate(data, replicates=100)


def test_default_controls_and_no_automatic_admission():
    assert module.calculate.__kwdefaults__ == {"replicates": 100000, "seed": 20260925}
    data = fixture()
    result = module.calculate(data)
    assert result["replicates"] == 100000 and result["seed"] == 20260925
    assert result["protocol_controls_match"] is True
    mean_family_f1 = np.mean([r["statistics_with_prior"]["F1"] for r in data["arms"][0]["families"]])
    assert abs(result["point_estimates"]["control"]["F1"] - mean_family_f1) > 1e-5
    assert result["publication_ready"] is False and result["uncertainty_admitted"] is False


@pytest.mark.parametrize("kwargs", [{"replicates": True}, {"replicates": 99}, {"seed": True}, {"seed": -1}])
def test_invalid_bootstrap_controls(kwargs):
    with pytest.raises(ValueError):
        module.calculate(fixture(), **kwargs)

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools.bootstrap_qfo_swiss_stages import aggregate, bootstrap, validated_values, METRICS


@pytest.fixture
def counts():
    return json.loads((Path(__file__).resolve().parents[2] /
                       "benchmark_tools/results/qfo_swiss_counts_20260917.json").read_text())


def test_aggregate_is_not_mean_f1():
    values = np.array([[0.9, 0.5], [0.5, 0.9]])
    assert aggregate(values.mean(axis=0))[0] == pytest.approx(0.7)
    assert aggregate(values).mean(axis=0)[0] != pytest.approx(0.7)


def test_shared_draws_zero_for_identical_stages(counts):
    for stage in counts["stages"][1:]:
        stage["families"] = copy.deepcopy(counts["stages"][0]["families"])
        stage["aggregate"] = counts["stages"][0]["aggregate"].copy()
    result = bootstrap(counts, replicates=100)
    for contrast in result["comparisons"]:
        for metric in contrast["metrics"].values():
            assert metric["difference"] == 0
            assert metric["paired_percentile_ci"] == [0, 0]
            assert metric["bonferroni_percentile_ci"] == [0, 0]
            assert metric["family_ties"] == 18


def test_direct_draw_reconstruction(counts):
    result = bootstrap(counts, replicates=100, seed=7)
    assert result == bootstrap(counts, replicates=100, seed=7)
    values = validated_values(counts)
    draws = np.random.Generator(np.random.PCG64(7)).multinomial(18, np.full(18, 1 / 18), size=100)
    diffs = []
    for draw in draws:
        indices = np.repeat(np.arange(18), draw)
        diffs.append(aggregate(values[2, indices].mean(axis=0)) - aggregate(values[0, indices].mean(axis=0)))
    for j, metric in enumerate(METRICS):
        row = result["comparisons"][0]["metrics"][metric]
        assert row["paired_percentile_ci"] == pytest.approx(np.quantile(np.asarray(diffs)[:, j], [0.025, 0.975]))
        assert row["bonferroni_percentile_ci"] == pytest.approx(np.quantile(np.asarray(diffs)[:, j], [0.05 / 24, 1 - 0.05 / 24]))
        assert row["family_wins"] + row["family_ties"] + row["family_losses"] == 18


@pytest.mark.parametrize("bad", ["stage", "family", "negative", "fractional", "duplicate", "overlap", "stat", "aggregate", "truth"])
def test_validation_failures(counts, bad):
    stage = counts["stages"][0]
    row = stage["families"][0]
    if bad == "stage":
        counts["stages"].reverse()
    elif bad == "family":
        stage["families"].reverse()
    elif bad in ("negative", "fractional"):
        row["counts_without_prior"]["TP"] = -1 if bad == "negative" else 0.5
    elif bad == "duplicate":
        row["represented_genes"].append(row["represented_genes"][0])
    elif bad == "overlap":
        stage["families"][1]["represented_genes"].append(row["represented_genes"][0])
    elif bad == "stat":
        row["statistics_with_prior"]["PPV"] = float("nan")
    elif bad == "aggregate":
        stage["aggregate"]["F1"] += 0.1
    else:
        counts["reference_relation_count"] += 1
    with pytest.raises(ValueError):
        bootstrap(counts, replicates=100)

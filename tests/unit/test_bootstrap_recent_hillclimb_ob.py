from copy import deepcopy

import numpy as np
import pytest

from benchmark_tools import bootstrap_recent_hillclimb_ob as reader
from benchmark_tools.bootstrap_orthobench import METRICS, statistics, weighted_records


def cell():
    records = [{"refog": f"RefOG{i:03}.txt", "genes": 4, "true_positive": 6,
                "false_positive": 0, "false_negative": 0, "exact": True, "splits": 1} for i in range(70)]
    score = {"refog_records": records, **dict.fromkeys(METRICS, 100.)}
    return {"dataset": "OrthoBench", "status": "independent_content_readback_verified",
            "coverage": {"genes": 251378, "species": 12}, "scores": {"OrthoBench: OrthoBench F1": 1.},
            "details": {"family_details": score}, "checked": []}


def panel():
    control = cell()
    first, second = deepcopy(control), deepcopy(control)
    first["matched_control"] = {"cell": deepcopy(control)}
    second["matched_control"] = deepcopy(control)
    return {"methods": [{"key": key, "scores": data["scores"], "cells": [data]}
                        for key, data in zip(reader.METHODS, (first, second))]}


def test_complete_family_records_reproduce_points_and_matched_controls():
    records, points, pins = reader.collect(panel())
    assert set(records) == {"matched_control", *reader.METHODS}
    assert all(np.array_equal(value, [100., 100., 100.]) for value in points.values())
    assert pins == []


@pytest.mark.parametrize("status", ["running", "failed", "not_run_failed_control_gate"])
def test_unverified_candidate_is_not_analyzed(status):
    summary = panel()
    summary["methods"][0]["cells"][0]["status"] = status
    with pytest.raises(ValueError, match="verified OrthoBench"):
        reader.collect(summary)


def test_equal_aggregate_control_scores_do_not_certify_family_identity():
    summary = panel()
    summary["methods"][1]["cells"][0]["matched_control"]["details"]["family_details"]["refog_records"][0]["splits"] = 2
    with pytest.raises(ValueError, match="different RefOG records"):
        reader.collect(summary)


@pytest.mark.parametrize("metric", METRICS)
def test_inconsistent_sufficient_statistics_rejected(metric):
    value = cell()
    value["details"]["family_details"][metric] -= 1
    with pytest.raises(ValueError, match="do not reproduce"):
        reader.sufficient_statistics(value)


def test_exactly_seventy_families_required():
    value = cell()
    value["details"]["family_details"]["refog_records"].pop()
    with pytest.raises(ValueError, match="70 RefOGs"):
        reader.sufficient_statistics(value)


def test_changed_gene_universe_or_control_status_rejected():
    value = cell()
    value["coverage"]["genes"] -= 1
    with pytest.raises(ValueError, match="complete OrthoBench universe"):
        reader.sufficient_statistics(value)
    value = cell()
    value["status"] = "failed"
    with pytest.raises(ValueError, match="complete OrthoBench universe"):
        reader.sufficient_statistics(value)


def test_admitted_method_and_cell_f1_must_agree():
    summary = panel()
    summary["methods"][0]["scores"] = {"OrthoBench: OrthoBench F1": .99}
    with pytest.raises(ValueError, match="verified OrthoBench"):
        reader.collect(summary)


def test_invalid_family_counts_never_generate_intervals():
    value = cell()
    value["details"]["family_details"]["refog_records"][0]["true_positive"] = -1
    with pytest.raises(ValueError):
        reader.sufficient_statistics(value)


def test_bootstrap_recomputes_weighted_statistics_not_mean_family_f1():
    value = cell()
    records = value["details"]["family_details"]["refog_records"]
    records[0].update(genes=20, true_positive=100, false_negative=90, false_positive=40)
    _, _, counts = weighted_records(records)
    aggregate = statistics(counts.sum(axis=0))
    value["details"]["family_details"].update(dict(zip(METRICS, aggregate)))
    value["scores"]["OrthoBench: OrthoBench F1"] = aggregate[0] / 100
    _, actual = reader.sufficient_statistics(value)
    assert np.array_equal(actual, aggregate)
    assert actual[0] != pytest.approx(statistics(counts)[:, 0].mean())

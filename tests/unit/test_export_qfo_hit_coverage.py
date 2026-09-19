import copy

import numpy as np
import pytest

from benchmark_tools.export_qfo_hit_coverage import validate
from benchmark_tools.summarize_sorted_search_hits import summarize, overlap


def fixture():
    arrays = {
        "hmm": ([0, 0, 1], [0, 1, 1]),
        "all_hits": ([0, 0, 1, 1], [0, 1, 0, 1]),
        "top100": ([0, 1], [0, 1]),
    }
    arrays = {k: (np.array(q), np.array(t), np.ones(len(q))) for k, (q, t) in arrays.items()}
    return {"status": "corrected_qfo_label_free_hit_comparison_complete",
            "accuracy_evaluated": False, "publication_ready": False,
            "species_labels": ["a", "b"],
            "searches": {k: summarize(*v, np.array([0, 1])) for k, v in arrays.items()},
            "overlaps": {a + "_vs_" + b: overlap(arrays[a], arrays[b], 2)
                         for a, b in (("hmm", "all_hits"), ("hmm", "top100"), ("all_hits", "top100"))}}


def test_valid_summary_is_not_mutated():
    report = fixture()
    original = copy.deepcopy(report)
    validate(report)
    assert report == original


@pytest.mark.parametrize("problem", ["accuracy", "status", "genes", "hits", "directions",
                                    "histogram", "overlap", "fraction", "coverage"])
def test_inconsistent_summary_rejected(problem):
    report = fixture()
    row = report["searches"]["hmm"]
    if problem == "accuracy":
        report["accuracy_evaluated"] = True
    elif problem == "status":
        report["status"] = "failed"
    elif problem == "genes":
        row["genes"] += 1
    elif problem == "hits":
        row["nonself_hits"] += 1
    elif problem == "directions":
        row["species_directions"].pop()
    elif problem == "histogram":
        row["normalized_score_histogram"]["counts"][0] += 1
    elif problem == "overlap":
        report["overlaps"]["hmm_vs_all_hits"]["all"]["intersection"] += 1
    elif problem == "fraction":
        report["overlaps"]["hmm_vs_all_hits"]["nonself"]["jaccard"] = float("nan")
    else:
        row["queries_without_hits"] = row["genes"] + 1
    with pytest.raises(ValueError):
        validate(report)

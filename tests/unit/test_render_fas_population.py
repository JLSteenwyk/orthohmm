from copy import deepcopy

import pytest

from benchmark_tools.render_fas_population import render, validate


def fixture():
    rows, methods = [], []
    for i in range(8):
        methods.append(dict(key=str(i), label="Method " + str(i), details=dict(FAS=dict(assessed_relations=10))))
        rows.append(dict(method=str(i), native_counts_match=True, saved_lookup_strata_and_values_match=True,
            database_historically_hash_bound=False, distinct_query_pairs=13, skipped_alias_pairs=1,
            precomputed=8, missing=2, unannotated=2, eligible_pairs=10, precomputed_score_sum=4.,
            native_logged_counts=dict(precomputed=8, missing=2, unannotated=2),
            hypothetical_full_mean_bounds=[.4, .6]))
    return dict(status="retained_fas_eligible_populations_recounted", methods=rows,
                uncertainty_admitted=False, benchmark_scores_changed=False, publication_ready=False), dict(methods=methods)


def test_render_scope_and_complete_order():
    report, manifest = fixture()
    text = render(report, manifest)
    assert "not\nconfidence intervals" in text
    assert text.count("| Method ") == 9
    assert text.index("| Method 0 |") < text.index("| Method 7 |")
    assert "0.400000 - 0.600000" in text


@pytest.mark.parametrize("field", ["uncertainty_admitted", "benchmark_scores_changed", "publication_ready"])
def test_rejects_broad_admission_flags(field):
    report, manifest = fixture()
    report[field] = True
    with pytest.raises(ValueError):
        validate(report, manifest)


@pytest.mark.parametrize("field", ["native_counts_match", "saved_lookup_strata_and_values_match"])
def test_rejects_unverified_recounts(field):
    report, manifest = fixture()
    report["methods"][0][field] = False
    with pytest.raises(ValueError):
        validate(report, manifest)


@pytest.mark.parametrize("change", ["partial", "duplicate", "reorder"])
def test_never_renders_selected_passing_subset(change):
    report, manifest = fixture()
    if change == "partial":
        report["methods"].pop()
    elif change == "duplicate":
        report["methods"][7] = deepcopy(report["methods"][0])
    else:
        report["methods"].reverse()
    with pytest.raises(ValueError):
        validate(report, manifest)


@pytest.mark.parametrize("field,value", [("eligible_pairs", 11), ("precomputed", 9),
    ("distinct_query_pairs", 12), ("precomputed_score_sum", float("nan")),
    ("precomputed_score_sum", 8.1), ("database_historically_hash_bound", True),
    ("missing", True), ("hypothetical_full_mean_bounds", [.4, .61])])
def test_rejects_inconsistent_arithmetic_or_provenance(field, value):
    report, manifest = fixture()
    report["methods"][0][field] = value
    with pytest.raises(ValueError):
        validate(report, manifest)

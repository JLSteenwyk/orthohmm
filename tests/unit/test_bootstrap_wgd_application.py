from copy import deepcopy

import pytest

from benchmark_tools.bootstrap_wgd_application import CONTRASTS, ENDPOINTS, intervals


def fixture():
    cohort = [{"orf_pair": [f"a{i}", f"b{i}"], "reference_eligible": True,
               "reference_pillar": "p" if i < 2 else "q",
               "available_members_by_species": {"Scerevisiae": 2, "Smikatae": 1}} for i in range(3)]
    methods = {}
    for name in {m for pair in CONTRASTS for m in pair}:
        methods[name] = [{**deepcopy(row), "coverage_denominator": 1,
                          **{e: int(name == "orthohmm_satellite_v2" and i < 2) for e in ENDPOINTS}}
                         for i, row in enumerate(cohort)]
    return cohort, methods


def test_pair_weighted_statistic_and_pillar_resampling():
    result = intervals(*fixture())
    assert result["population_pairs"] == 3 and result["population_pillars"] == 2
    assert len(result["comparisons"]) == result["multiplicity"] == 12
    row = result["comparisons"][0]
    assert row["difference_pp"] == pytest.approx(200 / 3)
    assert row["wins"] == 2 and row["ties"] == 1
    assert row["nominal95_pp"] == row["bonferroni12_pp"] == [0, 100]


def test_order_invariance_and_deterministic_shared_draws():
    cohort, methods = fixture()
    expected = intervals(cohort, methods)
    assert intervals(cohort[::-1], {k: v[::-1] for k, v in methods.items()}) == expected
    assert expected["comparisons"][0]["nominal95_pp"] == expected["comparisons"][1]["nominal95_pp"]


def test_failed_method_not_imputed_or_dropped_from_correction():
    cohort, methods = fixture()
    methods["sonicparanoid"] = None
    report = intervals(cohort, methods)
    assert report["multiplicity"] == 12
    assert sum(r["status"] == "unavailable_method" for r in report["comparisons"]) == 3


@pytest.mark.parametrize("mutation", ["missing_row", "missing_endpoint", "denominator", "fractional_split"])
def test_method_dependent_changes_rejected(mutation):
    cohort, methods = fixture()
    rows = methods["orthohmm_satellite_v2"]
    if mutation == "missing_row":
        rows.pop()
    elif mutation == "missing_endpoint":
        rows[0]["separation_rate"] = None
    elif mutation == "denominator":
        rows[0]["coverage_denominator"] = 2
    else:
        rows[0]["separation_rate"] = 0.5
    with pytest.raises(ValueError):
        intervals(cohort, methods)


def test_undefined_coverage_draws_are_not_silently_discarded():
    cohort, methods = fixture()
    cohort[2]["available_members_by_species"] = {"Scerevisiae": 2}
    for rows in methods.values():
        rows[2]["coverage_denominator"] = 0
        rows[2]["mean_non_scer_coverage"] = None
    result = intervals(cohort, methods)
    row = result["comparisons"][2]
    assert row["status"] == "undefined_replicates" and row["undefined_replicates"] > 0
    assert row["nominal95_pp"] is None and row["bonferroni12_pp"] is None

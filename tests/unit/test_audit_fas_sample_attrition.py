import pytest

from benchmark_tools.audit_fas_sample_attrition import analyze, one, scoring_task, render_table


LOG = """4 pairs precomputed, 2 missing (will compute); 3 no feature annotations
we will compute 2 new pairs and sample 4 precomputed pairs
FAS score[precomputed]: 0.750000 +- 0.100000 [N=4]
FAS score[missing]: 0.250000 +- 0.000000 [N=1]
"""


def test_missing_score_bounds():
    value = analyze(LOG, [0.75] * 4 + [0.25], 6)
    assert value["omitted_requested_new_scores"] == 1
    assert value["intended_sample_size"] == 6
    assert value["intended_sample_mean_bounds"] == [3.25 / 6, 4.25 / 6]
    assert value["saved_mean"] == 0.65


def test_no_missing_scores_collapses_bound():
    value = analyze(LOG.replace("[N=1]", "[N=2]"), [0.75] * 4 + [0.25] * 2, 6)
    assert value["omitted_requested_new_scores"] == 0
    assert value["intended_sample_mean_bounds"] == [value["saved_mean"]] * 2


@pytest.mark.parametrize("old,new", [
    ("4 pairs precomputed", "5 pairs precomputed"),
    ("compute 2 new", "compute 1 new"),
    ("sample 4 precomputed", "sample 3 precomputed"),
    ("[N=4]", "[N=3]"),
    ("[N=1]", "[N=0]"),
    ("0.750000", "0.700000"),
])
def test_inconsistent_log_rejected(old, new):
    with pytest.raises(ValueError):
        analyze(LOG.replace(old, new), [0.75] * 4 + [0.25], 6)


@pytest.mark.parametrize("values", [[0.75] * 4, [0.75] * 4 + [0.25, 0.25],
                                    [float("nan")] * 4 + [0.25], [1.1] * 4 + [0.25]])
def test_invalid_raw_scores_rejected(values):
    with pytest.raises(ValueError):
        analyze(LOG, values, 6)


def test_ambiguous_or_missing_log_record_rejected():
    for text in ("", LOG + LOG):
        with pytest.raises(ValueError):
            one(r"(\d+) pairs precomputed", text)


def test_scoring_task_not_downstream_shared_results():
    name = "FAS_example-one_raw.txt.gz"
    assert scoring_task('#!/bin/bash -ue\nfas_benchmark.py --participant "example_one"', name)
    assert not scoring_task("plot.py --participant example_one", name)
    assert not scoring_task("fas_benchmark.py --participant other", name)


@pytest.mark.parametrize("command", ["fas_benchmark.py --participant", "fas_benchmark.py",
    "fas_benchmark.py --participant --cpus 4", "fas_benchmark.py --participant a --participant b"])
def test_bad_scoring_command_rejected(command):
    with pytest.raises(ValueError):
        scoring_task(command, "FAS_a_raw.txt.gz")


def test_table_uses_computed_bounds_with_qualification():
    row = dict(method="example", **analyze(LOG, [0.75] * 4 + [0.25], 6))
    table = render_table(dict(methods=[row]))
    assert "not confidence intervals" in table
    assert "| example | 2 | 1 | 1 | 0.650000 | 0.541667 - 0.708333 |" in table

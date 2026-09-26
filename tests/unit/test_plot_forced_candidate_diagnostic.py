import copy

import pytest

from benchmark_tools.plot_forced_candidate_diagnostic import counts


def report():
    return dict(status="forced_candidate_scores_independently_recounted",
        candidate_only_interpretation_authorized=True, numerical_disagreements=[],
        directions=144, directed_pairs=81466, previously_scored=33098,
        decisions=dict(accepted=61975, scored_not_significant=19491),
        forced_transitions={"accepted:accepted": 31479,
            "scored_not_significant:scored_not_significant": 1619,
            "not_selected_by_prefilter:accepted": 30496,
            "not_selected_by_prefilter:scored_not_significant": 17872})


def test_valid_counts():
    assert counts(report()) == [[31479, 0], [0, 1619], [30496, 17872]]


@pytest.mark.parametrize("field,value", [
    ("numerical_disagreements", [1]), ("candidate_only_interpretation_authorized", False),
    ("previously_scored", 1), ("directions", 143), ("directed_pairs", 1),
    ("decisions", dict(accepted=1, scored_not_significant=19491)),
])
def test_reject_invalid_marginals(field, value):
    data = copy.deepcopy(report())
    data[field] = value
    with pytest.raises(ValueError):
        counts(data)

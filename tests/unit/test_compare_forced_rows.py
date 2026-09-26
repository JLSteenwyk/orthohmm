import pytest

from benchmark_tools.audit_ob_search_decisions import compare_forced_rows


def row(decision="accepted", score=2.0, evalue=1e-6):
    return dict(query="a", target="b", decision=decision, score=score, evalue=evalue)


def test_numeric_identity():
    result = compare_forced_rows([row()], [row()])
    assert result["previously_scored"] == 1
    assert result["numerical_disagreements"] == []
    assert result["transitions"] == {"accepted:accepted": 1}


def test_rescued_exclusion_not_numerical_disagreement():
    result = compare_forced_rows([row("not_selected_by_prefilter", None, None)], [row()])
    assert result["previously_scored"] == 0
    assert result["numerical_disagreements"] == []
    assert result["transitions"] == {"not_selected_by_prefilter:accepted": 1}


@pytest.mark.parametrize("change", [dict(score=2.1), dict(evalue=2e-6), dict(decision="scored_not_significant")])
def test_retains_numerical_or_decision_disagreement(change):
    result = compare_forced_rows([row()], [{**row(), **change}])
    assert len(result["numerical_disagreements"]) == 1


@pytest.mark.parametrize("prior,current", [
    ([row(), row()], [row()]), ([row()], []),
    ([row()], [row("not_selected_by_prefilter", None, None)]),
])
def test_rejects_incomplete_comparison(prior, current):
    with pytest.raises(ValueError):
        compare_forced_rows(prior, current)

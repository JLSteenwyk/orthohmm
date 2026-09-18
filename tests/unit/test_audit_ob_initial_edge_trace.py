import pytest

from benchmark_tools.audit_ob_initial_edge_trace import audit_rows


def fixture(score="2", threshold="2", decision="initial_edge"):
    source = dict(refog="F", left="a", right="b", forward_normalized_hit=score,
                  reverse_normalized_hit="NA", multipass_refined="True", root_hogs="False")
    traced = dict(refog="F", left="a", right="b", forward_score=score, reverse_score="NA",
                  left_threshold=threshold, right_threshold="inf", decision=decision,
                  multipass_refined="True", root_hogs="False")
    return source, traced


@pytest.mark.parametrize("score,threshold,decision", [
    ("2", "2", "initial_edge"), ("2", "3", "below_endpoint_threshold"),
    ("2", "inf", "no_finite_endpoint_threshold"), ("NA", "2", "no_direct_hit")])
def test_decisions(score, threshold, decision):
    source, traced = fixture(score, threshold, decision)
    result = audit_rows([source], [traced])
    assert result["summary"] == {decision + "/root_False": 1}
    assert result["pair_memberships"] == 1


@pytest.mark.parametrize("field,value", [("forward_score", "3"), ("forward_score", "nan"),
    ("decision", "below_endpoint_threshold"), ("root_hogs", "True"),
    ("root_hogs", "invalid"), ("left_threshold", "nan"), ("left_threshold", "-1")])
def test_corruption(field, value):
    source, traced = fixture()
    traced[field] = value
    with pytest.raises(ValueError):
        audit_rows([source], [traced])


def test_missing_and_duplicate_pairs():
    source, traced = fixture()
    for rows in ([], [traced, traced]):
        with pytest.raises(ValueError):
            audit_rows([source], rows)


def test_inconsistent_threshold():
    source, traced = fixture()
    second_source = dict(source, right="c")
    second_trace = dict(traced, right="c", left_threshold="1")
    with pytest.raises(ValueError, match="Inconsistent"):
        audit_rows([source, second_source], [traced, second_trace])

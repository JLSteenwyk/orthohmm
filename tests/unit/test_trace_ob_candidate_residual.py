from copy import deepcopy

import pytest

from benchmark_tools.trace_ob_candidate_residual import event_key, compare_events


def event(source=None, target=None):
    return dict(iteration=0, source_genes=source or ["a"], target_genes=target or ["b", "c"],
                source_cluster=1, target_cluster=2, support=1., margin=2., forward_average=3.,
                reverse_average=3., forward_normalized_support=1., reverse_normalized_support=1.)


def test_identity_ignores_labels_and_gene_order_not_direction_or_iteration():
    a, b = event(), event(target=["c", "b"])
    b["source_cluster"] = 10
    assert event_key(a) == event_key(b)
    assert event_key(a) != event_key(event(source=["b", "c"], target=["a"]))
    b["iteration"] = 1
    assert event_key(a) != event_key(b)


@pytest.mark.parametrize("field,value", [("source_genes", []), ("source_genes", ["a", "a"]),
    ("source_genes", ["b"]), ("target_genes", [None]), ("iteration", -1), ("iteration", True)])
def test_invalid_event_rejected(field, value):
    item = event()
    item[field] = value
    with pytest.raises(ValueError):
        event_key(item)


def test_order_difference_is_not_membership_difference():
    a, b = event(), event(source=["d"])
    result = compare_events([a, b], [b, a])
    assert result["common_events"] == 2
    assert result["left_only"] == result["right_only"] == []
    assert result["first_semantic_order_difference_zero_based"] == 0


def test_changed_events_and_numerical_values():
    a = event()
    b = deepcopy(a)
    b["support"] += 1e-12
    result = compare_events([a, event(source=["d"])], [b, event(source=["e"])])
    assert result["common_events"] == 1
    assert len(result["left_only"]) == len(result["right_only"]) == 1
    assert result["first_semantic_order_difference_zero_based"] == 1
    assert result["common_event_numerical_differences"]["support"]["changed_events"] == 1


def test_length_difference_and_infinite_margins():
    a = event()
    a["margin"] = float("inf")
    assert compare_events([a], [a])["common_event_numerical_differences"]["margin"]["max_absolute_difference"] == 0
    assert compare_events([a], [a, event(source=["d"])])["first_semantic_order_difference_zero_based"] == 1


def test_duplicate_semantic_event_rejected():
    with pytest.raises(ValueError, match="Duplicate"):
        compare_events([event(), event()], [event()])


def test_nan_rejected():
    a = event()
    a["support"] = float("nan")
    with pytest.raises(ValueError, match="NaN"):
        compare_events([a], [event()])

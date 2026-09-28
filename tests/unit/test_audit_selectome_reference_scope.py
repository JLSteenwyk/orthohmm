import pytest

from benchmark_tools.audit_selectome_reference_scope import summarize


def test_outside_relations_count_once_and_keep_truth_classes():
    truth = {("pooled", "a", "b"): True, ("pooled", "b", "c"): False,
             ("pooled", "c", "d"): True}
    result = summarize(truth, {"a": "in", "b": "in", "c": "out", "d": "out"}, {"in"})
    assert result["incident_proteins"] == 4
    assert result["incident_proteins_outside"] == 2
    assert result["relations"] == dict(both_inside=1, both_inside_ortholog=1,
        at_least_one_outside=2, at_least_one_outside_paralog=1, at_least_one_outside_ortholog=1)


@pytest.mark.parametrize("owners", [{"a": "in"}, {"a": "in", "b": ""},
    {"a": "in", "b": "in", "extra": "out"}])
def test_incomplete_extra_or_empty_assignment_rejected(owners):
    with pytest.raises(ValueError):
        summarize({("pooled", "a", "b"): True}, owners, {"in"})

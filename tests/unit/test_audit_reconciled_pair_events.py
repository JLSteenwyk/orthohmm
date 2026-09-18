import pytest

from benchmark_tools.audit_reconciled_pair_events import reconstruct


def rows():
    return [
        dict(node_id="a", parent_node_id="r", event="leaf", pair_event="leaf", genes="a", species="A"),
        dict(node_id="b", parent_node_id="r", event="leaf", pair_event="leaf", genes="b", species="B"),
        dict(node_id="r", parent_node_id="", event="speciation", pair_event="speciation", genes="a,b",
             species="A,B", mapping_conflict="false", species_overlap_count="0")]


def test_speciation_and_membership_filter():
    data = rows()
    groups = {"a": "g1", "b": "g2"}
    assert reconstruct(data, groups, False)[0] == {("a", "b")}
    pairs, summary = reconstruct(data, groups, True)
    assert pairs == set() and summary["group_filtered_pairs"] == 1


def test_mapping_conflict_is_uncertain_not_pair_duplication():
    data = rows()
    data[-1].update(event="duplication", pair_event="uncertain", mapping_conflict="true")
    assert reconstruct(data, {"a": "g", "b": "g"}, False)[0] == {("a", "b")}


def test_species_overlap_blocks_cross_child_pairs():
    data = rows()
    data[1]["species"] = "A"
    data[-1].update(event="duplication", pair_event="duplication", species="A", species_overlap_count="1")
    assert not reconstruct(data, {"a": "g", "b": "g"}, False)[0]


def test_duplication_blocks_cross_species_pairs_between_subtrees():
    data = rows()
    data[0]["parent_node_id"] = data[1]["parent_node_id"] = "left"
    data[-1].update(node_id="left", parent_node_id="root")
    data.extend([
        dict(node_id="c", parent_node_id="right", event="leaf", pair_event="leaf", genes="c", species="A"),
        dict(node_id="d", parent_node_id="right", event="leaf", pair_event="leaf", genes="d", species="C"),
        dict(node_id="right", parent_node_id="root", event="speciation", pair_event="speciation", genes="c,d",
             species="A,C", mapping_conflict="false", species_overlap_count="0"),
        dict(node_id="root", parent_node_id="", event="duplication", pair_event="duplication", genes="a,b,c,d",
             species="A,B,C", mapping_conflict="false", species_overlap_count="1")])
    assert reconstruct(data, dict.fromkeys("abcd", "g"), False)[0] == {("a", "b"), ("c", "d")}


@pytest.mark.parametrize("field,value", [("genes", "a,a"), ("genes", "a"), ("species", "A"),
    ("event", "duplication"), ("pair_event", "uncertain"), ("mapping_conflict", "unknown"),
    ("species_overlap_count", "1"), ("parent_node_id", "missing")])
def test_invalid_node_fields(field, value):
    data = rows()
    data[-1][field] = value
    with pytest.raises(ValueError):
        reconstruct(data, {"a": "g", "b": "g"}, False)


def test_incomplete_groups():
    with pytest.raises(ValueError):
        reconstruct(rows(), {"a": "g"}, False)


def test_duplicate_node():
    data = rows()
    with pytest.raises(ValueError):
        reconstruct(data + [data[0]], {"a": "g", "b": "g"}, False)


def test_disconnected_cycle():
    data = rows()
    data.extend([dict(node_id="x", parent_node_id="y"), dict(node_id="y", parent_node_id="x")])
    with pytest.raises(ValueError, match="Disconnected"):
        reconstruct(data, {"a": "g", "b": "g"}, False)

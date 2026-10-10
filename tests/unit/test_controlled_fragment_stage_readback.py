"""Focused tests for the independent retained-stage reader."""

from collections import defaultdict
import json
from pathlib import Path

import pytest

from benchmark_tools import readback_controlled_fragment_stages as reader


def node(name, parent, genes, species, mapping="S0001", event="leaf", pair="leaf", conflict="false", overlap=0, confidence="not_applicable", support=""):
    return dict(node_id=name, parent_node_id=parent, genes=",".join(genes), species=",".join(species),
                species_tree_node=mapping, event=event, pair_event=pair, mapping_conflict=conflict,
                species_overlap_count=str(overlap), event_confidence=confidence, branch_support=support)


def simple_tree():
    return {"a": node("a", "r", ["a"], ["x"], mapping="x"),
            "b": node("b", "r", ["b"], ["y"], mapping="y"),
            "r": node("r", "", ["a", "b"], ["x", "y"], event="speciation", pair="speciation", confidence="high")}


def test_pin_accepts_absolute_only_and_rejects_changed_bytes(tmp_path):
    path = tmp_path / "data"
    path.write_text("test")
    ref = reader.pin(path)
    assert reader.check(dict(absolute_path=ref["path"], bytes=ref["bytes"], sha256=ref["sha256"])) == path
    path.write_text("changed")
    with pytest.raises(ValueError, match="Changed input"):
        reader.check(ref)


def test_partition_requires_complete_disjoint_groups():
    assert reader.partition({"one": ["a"], "two": ["b"]}, {"a", "b"}) == {"a": "one", "b": "two"}
    for groups in ({"one": []}, {"one": ["a"]}, {"one": ["a", "a", "b"]}, {"one": ["a", "b", "c"]}):
        with pytest.raises(ValueError):
            reader.partition(groups, {"a", "b"})


def test_bfs_separates_indirect_connectivity_from_direct_edges(tmp_path):
    path = tmp_path / "edges"
    path.write_text("a\tb\t1\nb\tc\t2\n")
    edges, components = reader.graph(path, {"a", "b", "c", "d"})
    assert ("a", "c") not in edges
    assert components["a"] == components["c"] != components["d"]


@pytest.mark.parametrize("text", ["a\tb\tnan\n", "a\tb\t0\n", "a\tz\t1\n", "a\tb\t1\nb\ta\t1\n"])
def test_graph_rejects_invalid_evidence(tmp_path, text):
    path = tmp_path / "edges"
    path.write_text(text)
    with pytest.raises(ValueError):
        reader.graph(path, {"a", "b"})


def test_bottom_up_partition_and_pair_ancestor():
    nodes = simple_tree()
    groups, high = reader.tree_partition(nodes, {"a", "b"}, {"a": "x", "b": "y"})
    assert groups == [{"a", "b"}] and high == {("a", "b")}
    assert reader.ancestor(nodes, "a", "b")["node_id"] == "r"


def test_mapping_conflict_is_uncertain_not_positive_paralogy():
    nodes = simple_tree()
    nodes["a"]["species_tree_node"] = "S0001"
    nodes["r"].update(event="duplication", pair_event="uncertain", mapping_conflict="true", event_confidence="medium")
    groups, high = reader.tree_partition(nodes, {"a", "b"}, {"a": "x", "b": "y"})
    assert groups == [{"a", "b"}] and high == set()


@pytest.mark.parametrize("support,confidence", [("", "medium"), ("0.9", "high")])
def test_root_overlap_splits_with_support_dependent_confidence(support, confidence):
    nodes = simple_tree()
    nodes["b"].update(species="x", species_tree_node="x")
    nodes["r"].update(species="x", species_tree_node="S0000", event="duplication", pair_event="duplication",
                       species_overlap_count="1", event_confidence=confidence, branch_support=support)
    groups, high = reader.tree_partition(nodes, {"a", "b"}, {"a": "x", "b": "x"})
    assert {frozenset(g) for g in groups} == {frozenset({"a"}), frozenset({"b"})} and not high


def test_nonroot_overlap_does_not_force_root_partition():
    nodes = simple_tree()
    nodes["b"].update(species="x", species_tree_node="x")
    nodes["r"].update(species="x", event="duplication", pair_event="duplication", species_overlap_count="1", event_confidence="medium")
    groups, high = reader.tree_partition(nodes, {"a", "b"}, {"a": "x", "b": "x"})
    assert groups == [{"a", "b"}] and not high


def test_cyclic_nodes_are_rejected():
    nodes = simple_tree()
    nodes["a"]["parent_node_id"] = "a"
    with pytest.raises(ValueError):
        reader.tree_partition(nodes, {"a", "b"}, {"a": "x", "b": "y"})


@pytest.mark.parametrize("supported", [True, False])
def test_logged_constraints_require_high_confidence_cross_side_support(supported):
    event = dict(source_genes=["a"], target_genes=["b"])
    groups, details = reader.constrain([{"a", "b"}], {("a", "b")} if supported else set(), [(7, event)])
    assert len(groups) == (1 if supported else 2)
    assert details[0]["supported"] is supported and details[0]["event_index"] == 7
    assert details[0]["supporting_pair"] == (["a", "b"] if supported else None)


def test_native_pair_filter_activation_is_not_inferred_from_roothogs():
    context = dict(seeds={"a": "s", "b": "s"}, predictions={("a", "b")}, hits=defaultdict(list), edges=set(),
                   components={"a": 1, "b": 1}, candidates={"a": "f", "b": "f"}, roots={"a": "r1", "b": "r2"},
                   nodes={"f": simple_tree()}, before={"f": {"a": 1, "b": 2}}, detachments={"f": []}, active=False)
    row = reader.observation(reader.METHODS[1], context, ("a", "b"))
    assert row["predicted"] and not row["same_root_hog"] and row["observed_location"] == "event_rule_retention"
    context["active"] = True
    with pytest.raises(ValueError, match="Native event/root rule"):
        reader.observation(reader.METHODS[1], context, ("a", "b"))


def test_orthofinder_unavailable_fields_are_null_not_absent_hits():
    row = reader.observation(reader.METHODS[2], dict(seeds={"a": "s", "b": "s"}, predictions=set()), ("a", "b"))
    assert row["observed_location"] == "within_mcl_native_pair_exclusion"
    assert row["search_status"] == "unavailable_matched_adapter"
    assert all(row[k] is None for k in ("hit_forward", "hit_reverse", "graph_direct", "graph_connected", "pair_event"))


def test_readback_does_not_import_the_producer_or_reconstruction():
    text = Path(reader.__file__).read_text()
    assert "import trace_controlled_fragment_stages" not in text
    assert "import reconstruct_reconciliation_trace" not in text


def test_changed_report_is_rejected_before_native_reads(tmp_path):
    path = tmp_path / "report.json"
    path.write_text(json.dumps({}))
    with pytest.raises(ValueError, match="Changed stage report"):
        reader.verify(path, "0" * 64)

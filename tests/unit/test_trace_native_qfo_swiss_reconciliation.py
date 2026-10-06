"""Topology and candidate-identity checks for observed native event localization."""

import copy
import csv
import hashlib

import pytest

from benchmark_tools import trace_native_qfo_swiss_reconciliation as trace


@pytest.fixture
def tree():
    nodes = {}
    for name, parent, species in (("a", "G0", "A"), ("b", "G0", "B"), ("c", "G1", "A"), ("d", "G1", "B")):
        nodes[name] = dict(source_family="Family0000000", node_id=name, parent_node_id=parent,
            event="leaf", species_tree_node=species, species={species}, genes={name}, pair_event="leaf",
            event_confidence="not_applicable", species_overlap_count=0, mapping_conflict="false", branch_support=None)
    for name, parent, genes, event in (("G0", "G2", {"a", "b"}, "speciation"),
                                     ("G1", "G2", {"c", "d"}, "speciation"),
                                     ("G2", "", {"a", "b", "c", "d"}, "duplication")):
        nodes[name] = dict(source_family="Family0000000", node_id=name, parent_node_id=parent,
            event=event, species_tree_node="S0", species={"A", "B"}, genes=genes, pair_event=event,
            event_confidence="high", species_overlap_count=2 if event == "duplication" else 0,
            mapping_conflict="true" if event == "duplication" else "false", branch_support=.99)
    return nodes


def test_descendant_validation_and_actual_lca(tree):
    trace.validate_nodes(tree, set("abcd"))
    assert trace.pair_lca(tree, "a", "c")["node_id"] == "G2"
    assert trace.pair_lca(tree, "a", "b")["pair_event"] == "speciation"


@pytest.mark.parametrize("fault", ("missing_parent", "self_parent", "second_root", "missing_leaf",
    "descendant_genes", "descendant_species", "overlap", "mapping", "pair_event", "event", "support_nan",
    "support_range", "bad_flag", "bad_confidence", "disconnected_cycle", "unary"))
def test_invalid_annotation_refuses(tree, fault):
    nodes = copy.deepcopy(tree)
    if fault == "missing_parent":
        nodes["a"]["parent_node_id"] = "absent"
    elif fault == "self_parent":
        nodes["a"]["parent_node_id"] = "a"
    elif fault == "second_root":
        nodes["a"]["parent_node_id"] = ""
    elif fault == "missing_leaf":
        nodes.pop("a")
    elif fault == "descendant_genes":
        nodes["G2"]["genes"] = set("abc")
    elif fault == "descendant_species":
        nodes["G2"]["species"] = {"A"}
    elif fault == "overlap":
        nodes["G2"]["species_overlap_count"] = 1
    elif fault == "mapping":
        nodes["G2"]["mapping_conflict"] = "false"
    elif fault == "pair_event":
        nodes["G2"]["pair_event"] = "speciation"
    elif fault == "event":
        nodes["G2"]["event"] = "speciation"
    elif fault == "support_nan":
        nodes["G2"]["branch_support"] = float("nan")
    elif fault == "support_range":
        nodes["G2"]["branch_support"] = 1.1
    elif fault == "bad_flag":
        nodes["G2"]["mapping_conflict"] = "maybe"
    elif fault == "bad_confidence":
        nodes["G2"]["event_confidence"] = "certain"
    elif fault == "disconnected_cycle":
        nodes["x"] = dict(nodes["a"], node_id="x", parent_node_id="y")
        nodes["y"] = dict(nodes["a"], node_id="y", parent_node_id="x")
    elif fault == "unary":
        nodes["b"]["parent_node_id"] = "G2"
    with pytest.raises(ValueError):
        trace.validate_nodes(nodes, set("abcd"))


def test_uncertain_mapping_conflict_is_not_paralog_exclusion(tree):
    nodes = copy.deepcopy(tree)
    nodes["c"]["species"] = {"C"}
    nodes["c"]["species_tree_node"] = "C"
    nodes["d"]["species"] = {"D"}
    nodes["d"]["species_tree_node"] = "D"
    nodes["G1"]["species"] = {"C", "D"}
    root = nodes["G2"]
    root.update(species=set("ABCD"), species_overlap_count=0, pair_event="uncertain", event_confidence="medium")
    trace.validate_nodes(nodes, set("abcd"))
    assert trace.pair_lca(nodes, "a", "c")["pair_event"] == "uncertain"


def test_streamed_root_reconstruction_and_target_map(tmp_path):
    path = tmp_path / "roots.tsv"
    path.write_text("root_hog\tsource_family\tgenes\n"
        "RootHOG0000000\tFamily0000000\ta,b\nRootHOG0000001\tFamily0000000\tc,d\n"
        "RootHOG0000002\tFamily0000001\te\n")
    digest = hashlib.sha256(b"a b c d\ne\n").hexdigest()
    selected, families, summary = trace.source_groups(path, {"a", "c"}, digest)
    assert families == {"Family0000000": set("abcd")}
    assert selected["a"]["root_hog"] != selected["c"]["root_hog"]
    assert summary == dict(source_families=2, root_hogs=3, input_genes=5, reconstructed_candidate_sha256=digest)
    with pytest.raises(ValueError, match="candidate groups"):
        trace.source_groups(path, {"a", "c"}, "0" * 64)
    with pytest.raises(ValueError, match="endpoint missing"):
        trace.source_groups(path, {"unknown"}, digest)


@pytest.mark.parametrize("text", (
    "wrong\n", "root_hog\tsource_family\tgenes\n",
    "root_hog\tsource_family\tgenes\nRootHOG0000001\tFamily0000000\ta\n",
    "root_hog\tsource_family\tgenes\nRootHOG0000000\tFamily0000001\ta\n",
    "root_hog\tsource_family\tgenes\nRootHOG0000000\tFamily0000000\ta,a\n",
    "root_hog\tsource_family\tgenes\nRootHOG0000000\tFamily0000000\tsp|A|x,tr|A|y\n"))
def test_bad_roots_refuse(tmp_path, text):
    path = tmp_path / "roots.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        trace.source_groups(path, {"A", "a"}, "0" * 64)


def test_structured_node_readback(tmp_path, tree):
    path = tmp_path / "nodes.tsv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=trace.NODE_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for node in tree.values():
            writer.writerow(dict(node, genes=",".join(sorted(node["genes"])), species=",".join(sorted(node["species"]))))
    nodes, count = trace.read_nodes(path, {"Family0000000": set("abcd")})
    assert count == 7 and nodes["Family0000000"] == tree


def test_native_pair_membership_stream(tmp_path):
    path = tmp_path / "pairs.tsv"
    path.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\na\tA\tb\tB\n")
    assert trace.selected_predictions(path, {("a", "b"), ("a", "c")}) == ({("a", "b")}, 1)
    path.write_text("wrong\n")
    with pytest.raises(ValueError):
        trace.selected_predictions(path, set())


def test_missing_lca_endpoint_refuses(tree):
    with pytest.raises(ValueError, match="endpoints missing"):
        trace.pair_lca(tree, "a", "absent")

import copy

import pytest

from benchmark_tools.reconstruct_reconciliation_trace import apply_logged_constraints, reconstruct_bypass, reconstruct_nodes


def node(key, parent, genes, species, event="leaf", overlap=0, mapped="S0000", conflict="false", confidence="not_applicable"):
    return {"node_id": key, "parent_node_id": parent, "genes": ",".join(genes), "species": ",".join(species),
            "event": event, "pair_event": "uncertain" if conflict == "true" and not overlap else event,
            "event_confidence": confidence, "mapping_conflict": conflict, "species_overlap_count": str(overlap), "species_tree_node": mapped}


def tree():
    owners = {"a": "x", "b": "x", "c": "y", "d": "z"}
    rows = [node("a", "D", ["a"], ["x"]), node("b", "D", ["b"], ["x"]),
            node("D", "R", ["a", "b"], ["x"], "duplication", 1, confidence="medium"),
            node("c", "S", ["c"], ["y"]), node("d", "S", ["d"], ["z"]),
            node("S", "R", ["c", "d"], ["y", "z"], "speciation", confidence="high"),
            node("R", "", ["a", "b", "c", "d"], ["x", "y", "z"], "speciation", confidence="high")]
    return rows, owners


def test_root_duplication_and_propagated_split_are_distinguished():
    rows, owners = tree()
    result = reconstruct_nodes(rows, set(owners), owners)
    assert {frozenset(g) for g in result["root_groups"]} == {frozenset("a"), frozenset("b"), frozenset("cd")}
    assert result["root_duplication_nodes"] == ["D"]
    assert result["propagated_split_nodes"] == ["R"]
    assert ("c", "d") in result["high_confidence_pairs"]
    assert ("a", "b") not in result["high_confidence_pairs"]


def test_duplication_below_species_root_does_not_split_root_groups():
    rows, owners = tree()
    rows[2]["species_tree_node"] = "S0001"
    assert reconstruct_nodes(rows, set(owners), owners)["root_groups"] == [set(owners)]


@pytest.mark.parametrize("change", ["duplicate", "missing_parent", "overlap", "species", "parent_genes", "event", "confidence", "disconnected"])
def test_invalid_node_evidence_rejected(change):
    rows, owners = tree()
    if change == "duplicate":
        rows.append(copy.deepcopy(rows[0]))
    elif change == "missing_parent":
        rows[0]["parent_node_id"] = "missing"
    elif change == "overlap":
        rows[2]["species_overlap_count"] = "0"
    elif change == "species":
        rows[0]["species"] = "y"
    elif change == "parent_genes":
        rows[2]["genes"] = "a"
    elif change == "event":
        rows[2]["pair_event"] = "speciation"
    elif change == "confidence":
        rows[-1]["event_confidence"] = "medium"
    elif change == "disconnected":
        rows += [node("orphan", "orphan", ["a"], ["x"])]
    with pytest.raises(ValueError):
        reconstruct_nodes(rows, set(owners), owners)


def test_bypass_rejects_missing_ambiguous_tree():
    with pytest.raises(ValueError, match="Ambiguous"):
        reconstruct_bypass(set("abc"), {"a": "x", "b": "x", "c": "y"})
    result = reconstruct_bypass(set("ab"), {"a": "x", "b": "y"})
    assert result["high_confidence_pairs"] == {("a", "b")}


def test_constraints_detach_only_unsupported_sources():
    reconstruction = {"root_groups": [set("abcd")], "high_confidence_pairs": {("a", "b")}}
    constraints = [(5, {"source_genes": ["a"], "target_genes": ["b"]}),
                   (6, {"source_genes": ["c"], "target_genes": ["d"]})]
    groups, evidence = apply_logged_constraints(reconstruction, constraints, set("abcd"))
    assert {frozenset(g) for g in groups} == {frozenset("abd"), frozenset("c")}
    assert groups == [set("abd"), {"c"}]
    assert evidence[0]["supporting_pair"] == ("a", "b")
    assert evidence[1]["supported"] is False
    with pytest.raises(ValueError, match="multiple detached"):
        apply_logged_constraints(reconstruction, [constraints[1], constraints[1]], set("abcd"))


@pytest.mark.parametrize("newick", ["((a,c),(b,d));", "(((a,b),c),d);", "(a,b,c,d);"])
def test_reconstruction_matches_native_reconciler(tmp_path, newick):
    from dataclasses import asdict
    from orthohmm.phylogeny import parse_gene_tree, parse_species_tree, reconcile_gene_tree
    path = tmp_path / "species.nwk"
    path.write_text("[&R] (X,Y);")
    owners = {"a": "X", "b": "X", "c": "Y", "d": "Y"}
    result = reconcile_gene_tree(parse_gene_tree(newick), parse_species_tree(path, ["X.faa", "Y.faa"]).tree,
                                 owners, root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")
    rows = []
    for entry in result.nodes:
        row = asdict(entry)
        row.update(parent_node_id=row["parent_node_id"] or "", genes=",".join(row["genes"]),
                   species=",".join(row["species"]), mapping_conflict=str(row["mapping_conflict"]).lower(),
                   species_overlap_count=str(row["species_overlap_count"]))
        rows.append(row)
    reconstructed = reconstruct_nodes(rows, set(owners), owners)
    assert {frozenset(g) for g in reconstructed["root_groups"]} == {frozenset(g) for g in result.root_groups}
    assert reconstructed["high_confidence_pairs"] == {(a, b) for a, b, level in result.ortholog_pair_confidence if level == "high"}

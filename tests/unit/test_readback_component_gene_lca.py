import pytest

from benchmark_tools.readback_component_gene_lca import lca_from_parent_table, raw_support, reconstructed_feature, split_key


def row(identifier, parent, genes, event, overlap=0, conflict=False, pair_event=None):
    return {"node_id": identifier, "parent_node_id": parent, "genes": ",".join(genes),
            "event": event, "species_overlap_count": str(overlap), "mapping_conflict": str(conflict).lower(),
            "pair_event": pair_event or event}


def nodes():
    return {name: value for name, value in (
        ("a",row("a","ab",["a"],"leaf")), ("b",row("b","ab",["b"],"leaf")),
        ("c",row("c","root",["c"],"leaf")), ("d",row("d","root",["d"],"leaf")),
        ("ab",row("ab","root",["a","b"],"speciation")),
        ("root",row("root","",list("abcd"),"speciation")))}


def test_raw_support_independently_restores_split_identity_and_percentage_scale():
    universe, values = raw_support("((g00000000:1,g00000001:1)72:1,g00000002:1,g00000003:1);", tuple("abcd"))
    assert values == {split_key("ab",universe): .72}
    assert reconstructed_feature(("a","b"), nodes(), universe, values) == ("reconciled", "speciation", "speciation", False, "0", "high", False, "0_7_to_0_9")
    assert reconstructed_feature(("a","c"), nodes(), universe, values)[-2:] == (True,"unavailable")


def test_raw_support_rejects_conflicting_split_values():
    with pytest.raises(ValueError, match="Conflicting"):
        raw_support("((g00000000:1,g00000001:1)0.2:1,(g00000002:1,g00000003:1)0.8:1);", tuple("abcd"))


def test_parent_lca_is_independent_of_child_order_and_rejects_disconnected_or_cyclic_records():
    table = nodes()
    assert lca_from_parent_table(("b","a"),table)["node_id"] == "ab"
    assert lca_from_parent_table(("d","a"),table)["node_id"] == "root"
    table["ab"]["parent_node_id"] = "ab"
    with pytest.raises(ValueError, match="ancestry"):
        lca_from_parent_table(("a","b"),table)


def test_corrected_support_can_change_duplication_event_confidence_without_pair_acceptance():
    table = nodes()
    table["ab"] = row("ab","root",["a","b"],"duplication",overlap=1,pair_event="duplication")
    universe = frozenset("abcd")
    for support, label in ((.89,"medium"),(.9,"high")):
        feature = reconstructed_feature(("a","b"), table, universe, {split_key("ab",universe):support})
        assert feature[2] == "duplication" and feature[5] == label


def test_mapping_only_conflict_remains_uncertain_and_unavailable_support_is_not_invented():
    table = nodes()
    table["ab"] = row("ab","root",["a","b"],"duplication",conflict=True,pair_event="uncertain")
    feature = reconstructed_feature(("a","b"), table, frozenset("abcd"), {})
    assert feature[2] == "uncertain" and feature[5] == "medium" and feature[-1] == "unavailable"


@pytest.mark.parametrize("genes", [(),tuple("abcd"),tuple("abcde")])
def test_split_key_rejects_missing_or_nonproper_branches(genes):
    with pytest.raises(ValueError, match="proper"):
        split_key(genes, frozenset("abcd"))

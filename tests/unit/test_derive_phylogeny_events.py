import dendropy
import pytest

from benchmark_tools.derive_phylogeny_events import derive, constrain


def tree(text):
    return dendropy.Tree.get(data=text, schema="newick", preserve_underscores=True, rooting="force-rooted")


def test_speciation_oracle():
    rows, groups, pairs = derive(tree("((a,b),c);"), tree("((A,B),C);"),
                                 {"a": "A", "b": "B", "c": "C"}, "F")
    assert pairs == {("a", "b"): "high", ("a", "c"): "high", ("b", "c"): "high"}
    assert groups == [frozenset("abc")]
    assert rows[-1]["event"] == rows[-1]["pair_event"] == "speciation"
    assert rows[-1]["parent_node_id"] == ""


def test_mapping_conflict_is_not_positive_paralogy():
    rows, groups, pairs = derive(tree("((a,c),b);"), tree("((A,B),C);"),
                                 {"a": "A", "b": "B", "c": "C"}, "F")
    assert rows[-1]["event"] == "duplication"
    assert rows[-1]["pair_event"] == "uncertain"
    assert pairs == {("a", "c"): "high", ("a", "b"): "medium", ("b", "c"): "medium"}
    assert groups == [frozenset("abc")]


@pytest.mark.parametrize("support,confidence", [("", "medium"), ("95", "high"),
                                               ("0.91", "high"), ("nan", "medium")])
def test_root_overlap_splits_lineages(support, confidence):
    rows, groups, pairs = derive(tree(f"((a1,b1),(a2,c1)){support};"), tree("((A,B),C);"),
        {"a1": "A", "b1": "B", "a2": "A", "c1": "C"}, "F")
    assert set(groups) == {frozenset(["a1", "b1"]), frozenset(["a2", "c1"])}
    assert pairs == {("a1", "b1"): "high", ("a2", "c1"): "high"}
    assert rows[-1]["pair_event"] == "duplication"
    assert rows[-1]["event_confidence"] == confidence


def test_subroot_overlap_does_not_split_root_group():
    _, groups, pairs = derive(tree("((a1,a2),b);"), tree("(A,B);"),
                              {"a1": "A", "a2": "A", "b": "B"}, "F")
    assert groups == [frozenset(["a1", "a2", "b"])]
    assert pairs == {("a1", "b"): "high", ("a2", "b"): "high"}


def test_satellite_detachment_and_global_filter():
    groups = [frozenset("abc")]
    pairs = {("a", "b"): "medium", ("a", "c"): "high", ("b", "c"): "medium"}
    refined, retained, counts = constrain(groups, pairs,
        [dict(source_genes=["b"], target_genes=["a", "c"])], True)
    assert set(refined) == {frozenset("ac"), frozenset("b")}
    assert retained == {("a", "c"): "high"}
    assert counts["ortholog_pairs_removed"] == 2
    assert counts["detached_genes"] == 1
    assert constrain(refined, pairs, [], False)[1] == pairs
    assert constrain(refined, pairs, [], True)[1] == retained


def test_supported_constraint():
    groups, pairs, counts = constrain([frozenset("ab")], {("a", "b"): "high"},
                                     [dict(source_genes=["a"], target_genes=["b"])], True)
    assert groups == [frozenset("ab")]
    assert counts["supported_constraints"] == 1
    assert counts["detached_genes"] == 0


@pytest.mark.parametrize("source,target", [([], ["a"]), (["a"], ["a"]), (["x"], ["b"])])
def test_bad_constraints(source, target):
    with pytest.raises(ValueError):
        constrain([frozenset("ab")], {}, [dict(source_genes=source, target_genes=target)], True)

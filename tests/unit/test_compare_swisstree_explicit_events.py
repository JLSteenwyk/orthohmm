import pytest
from Bio.Phylo.PhyloXML import Clade, Events, Phylogeny

from benchmark_tools.compare_swisstree_explicit_events import compare_tree, event_label


@pytest.mark.parametrize("events,expected", [(None, "unknown"), (Events(duplications=0), "unknown"),
    (Events(duplications=1), "duplication"), (Events(speciations=1), "speciation"),
    (Events(duplications=1, speciations=1), "ambiguous"), (Events(type="speciation_or_duplication"), "unknown")])
def test_no_implicit_speciation(events, expected):
    assert event_label(Clade(events=events)) == expected


def test_pair_diagnostic_retains_unrooted_and_unknown():
    left = Clade(clades=[Clade(name="A_HUMAN_a"), Clade(name="B_MOUSE_b")])
    tree = Phylogeny(root=Clade(events=Events(duplications=1), clades=[left, Clade(name="C_RAT_c")]), rooted=False)
    result = compare_tree(tree, {"a", "b", "c", "d"}, {("a", "b"): True, ("a", "c"): True,
                                                        ("b", "c"): False, ("c", "d"): False})
    assert result["pair_counts"] == dict(unknown_reference_ortholog=1, duplication_reference_ortholog=1,
                                       duplication_reference_nonortholog=1, unmapped_pairs=1)
    assert result["rooted_attribute"] is False
    assert len(result["disagreements"]) == 1
    assert result["missing_mapping_candidates"] == ["d"]
    assert result["ancestral_labels_admitted"] is False


def test_ambiguous_suffix_not_selected():
    tree = Phylogeny(root=Clade(clades=[Clade(name="A_x"), Clade(name="B_x"), Clade(name="C_y")]))
    result = compare_tree(tree, {"x", "y"}, {("x", "y"): True})
    assert result["ambiguous_mapping_candidates"] == {"x": ["A_x", "B_x"]}
    assert result["pair_counts"] == {"unmapped_pairs": 1}


def test_duplicate_leaves_rejected():
    tree = Phylogeny(root=Clade(clades=[Clade(name="x"), Clade(name="x")]))
    with pytest.raises(ValueError):
        compare_tree(tree, {"x", "y"}, {("x", "y"): True})


def test_same_leaf_for_multiple_accessions_rejected():
    tree = Phylogeny(root=Clade(clades=[Clade(name="A_x"), Clade(name="y")]))
    with pytest.raises(ValueError):
        compare_tree(tree, {"A_x", "x", "y"}, {("x", "y"): True})

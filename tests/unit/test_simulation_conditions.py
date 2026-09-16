from io import StringIO

from Bio import Phylo
import pytest

from benchmark_tools.simulation_conditions import group_pairs, project_truth, ranked_names, score_pairs, select_removals


def tree(text="(((A,B),(C,D)),((E,F),(G,H))); "):
    return Phylo.read(StringIO(text), "newick")


def test_missingness_is_exact_order_independent_and_seeded():
    genes = {f"g{i}": "A" if i % 2 else "B" for i in range(23)}
    first = select_removals(genes, tree(), 20261001, "missing20")
    assert len(first["removed_genes"]) == 4
    assert first["removed_species"] == []
    assert first == select_removals(dict(reversed(list(genes.items()))), tree(), 20261001, "missing20")
    assert first["removed_genes"] != select_removals(genes, tree(), 20261002, "missing20")["removed_genes"]


def test_clade_rule_and_count_control():
    genes = {s + "_1": s for s in "ABCDEFGH"}
    uneven = select_removals(genes, tree(), 3, "uneven_taxa")
    assert uneven["removed_species"] == ["B", "C", "D"]
    control = select_removals(genes, tree(), 3, "taxon_count_control")
    assert len(control["removed_species"]) == 3
    assert len(control["retained_species"]) == 5
    # Sibling ordering must not choose a different clade.
    assert uneven == select_removals(genes, tree("(((H,G),(F,E)),((D,C),(B,A)));"), 3, "uneven_taxa")


def test_no_eligible_clade_is_not_silently_substituted():
    result = select_removals({s: s for s in "ABCDEFGH"}, tree("(A,B,C,D,E,F,G,H);"), 1, "uneven_taxa")
    assert result["status"] == "inapplicable"


def test_projection_only_removes_incident_true_pairs():
    genes = {"a": "A", "b": "B", "c": "C"}
    retained, pairs = project_truth([("a", "b"), ("b", "c"), ("a", "c")], genes,
                                   {"status": "selected", "removed_genes": ["b"]})
    assert retained == {"a": "A", "c": "C"}
    assert pairs == [("a", "c")]


def test_unrelated_family_false_positives_are_counted():
    genes = {"F1a": "A", "F1b": "B", "F2a": "A", "F2b": "B"}
    truth = [("F1a", "F1b"), ("F2a", "F2b")]
    result = score_pairs(group_pairs([list(genes)], genes), truth, genes)
    assert (result["tp"], result["fp"], result["fn"]) == (2, 2, 0)
    assert result["f1"] == pytest.approx(2 / 3)
    assert result["pair_endpoint_coverage"] == 1


def test_pair_orientation_deduplicated_and_missing_predictions_not_imputed():
    genes = {"a": "A", "b": "B", "c": "C"}
    result = score_pairs([("b", "a"), ("a", "b")], [("a", "b"), ("b", "c")], genes)
    assert (result["tp"], result["fp"], result["fn"]) == (1, 0, 1)
    assert result["duplicate_prediction_rows"] == 1
    assert result["pair_endpoint_coverage"] == 2 / 3
    empty = score_pairs([], [], genes)
    assert empty["undefined_ratios"] == ["f1", "precision", "recall"]


@pytest.mark.parametrize("pairs", [[("a", "missing")], [("a", "a")], [("a", "aa")], ["ab"]])
def test_invalid_pairs_rejected(pairs):
    with pytest.raises(ValueError):
        score_pairs(pairs, [], {"a": "A", "aa": "A", "b": "B"})


def test_duplicate_group_membership_rejected():
    with pytest.raises(ValueError):
        list(group_pairs([["a", "b"], ["a"]], {"a": "A", "b": "B"}))
    with pytest.raises(ValueError):
        ranked_names(["a", "a"], 1, "missing20")

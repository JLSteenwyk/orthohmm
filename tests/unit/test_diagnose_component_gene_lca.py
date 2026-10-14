from dataclasses import replace
import pytest

from benchmark_tools.diagnose_component_gene_lca import category_and_score, node_index, pair_features, support_bin, verify_totals
from orthohmm.phylogeny import parse_gene_tree, reconcile_gene_tree


@pytest.mark.parametrize("value,expected", [(None,"unavailable"), (0.,"below_0_5"), (.4999,"below_0_5"),
    (.5,"0_5_to_0_7"), (.7,"0_7_to_0_9"), (.9,"at_least_0_9"), (1.,"at_least_0_9")])
def test_support_bin_boundaries_are_prespecified(value, expected):
    assert support_bin(value) == expected


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1., 1.01])
def test_support_bin_rejects_non_normalized_support(value):
    with pytest.raises(ValueError, match="normalized"):
        support_bin(value)


def reconciled():
    return reconcile_gene_tree(parse_gene_tree("((a:1,b:1)0.72:1,c:2)0.99;"),
        parse_gene_tree("((a:1,b:1):1,c:2);"), {gene:gene for gene in "abc"},
        root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")


def test_pair_lca_features_use_measured_internal_support_but_not_root_incoming_support():
    indexed, ancestry = node_index(reconciled().nodes)
    assert pair_features(("a","b"), indexed, ancestry) == ("speciation", "speciation", False, "0", "high", False, "0_7_to_0_9")
    assert pair_features(("a","c"), indexed, ancestry) == ("speciation", "speciation", False, "0", "high", True, "unavailable")


@pytest.mark.parametrize("pair", [("a","a"), ("a","missing"), ("a","G00000")])
def test_pair_features_rejects_non_leaf_or_collapsed_pairs(pair):
    indexed, ancestry = node_index(reconciled().nodes)
    with pytest.raises(ValueError, match="distinct"):
        pair_features(pair, indexed, ancestry)


def test_node_index_rejects_missing_parent_cycle_or_multiple_roots():
    nodes = reconciled().nodes
    for altered in (nodes + (nodes[0],), (replace(nodes[0], parent_node_id="missing"), *nodes[1:]),
        (replace(nodes[0], parent_node_id=nodes[0].node_id), *nodes[1:]),
        (replace(nodes[0], parent_node_id=None), *nodes[1:])):
        with pytest.raises(ValueError):
            node_index(altered)


def test_pair_features_preserves_positive_paralogy_mapping_only_conflict():
    result = reconcile_gene_tree(parse_gene_tree("((a:1,c:1)0.95:1,b:2);"),
        parse_gene_tree("((a:1,b:1):1,c:2);"), {gene:gene for gene in "abc"},
        root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")
    indexed, ancestry = node_index(result.nodes)
    assert pair_features(("a","b"), indexed, ancestry)[:5] == ("duplication", "uncertain", True, "0", "medium")


def test_categories_reconstruct_added_removed_shared_scores_and_reject_changed_shared_values():
    a, b, c = ("a","b"), ("a","c"), ("b","c")
    left, right = {a:1000000,b:250000}, {a:1000000,c:500000}
    assert category_and_score(a,left,right) == ("shared",1000000)
    assert category_and_score(b,left,right) == ("removed",250000)
    assert category_and_score(c,left,right) == ("added",500000)
    with pytest.raises(ValueError, match="Shared"):
        category_and_score(a,left,{a:999999})


def test_bucket_verifier_requires_original_counts_and_integer_sums():
    rows = [{"metric": metric, "category": category, "scored_relations": 1, "score_sum_millionths": score}
            for metric in ("GO", "EC") for category, score in (("shared", 1000000), ("added", 500000), ("removed", 250000))]
    counts = {metric + ":" + cell: {"pairs":2, "sum_millionths":total}
              for metric in ("GO", "EC") for cell, total in (("control",1250000), ("components500",1500000))}
    verify_totals(rows, counts)
    for key in ("scored_relations", "score_sum_millionths"):
        altered = [dict(row) for row in rows]
        altered[0][key] += 1
        with pytest.raises(ValueError, match="counts/sums"):
            verify_totals(altered, counts)

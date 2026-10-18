from copy import deepcopy
from dataclasses import replace
import json

import numpy as np
import pytest

from orthohmm.local_rescue import assess_joint_outcome, nominate_local_contexts


def fixture(sizes=(3, 3, 3, 3), links=((0, 1, 9.), (0, 2, 4.), (2, 3, 6.))):
    names, species, clusters = [], [], []
    for family, size in enumerate(sizes):
        group = []
        for index in range(size):
            group.append(len(names))
            names.append(f"family{family}_gene{index}")
            species.append(index % 3)
        clusters.append(group)
    q, t, scores = [], [], []
    for a, b, value in links:
        aa, bb = clusters[a], clusters[b]
        for left, right in ((aa[0], bb[1]), (aa[1], bb[2]), (bb[1], aa[0]), (bb[2], aa[1])):
            q.append(left)
            t.append(right)
            scores.append(value)
    return (clusters, names, np.asarray(species, dtype=np.int32), np.asarray(q, dtype=np.int32),
            np.asarray(t, dtype=np.int32), np.asarray(scores, dtype=np.float64))


def semantic(result):
    value = deepcopy(result)
    value["parameters"].pop("hit_chunk_size")
    return value


def test_unique_mutual_best_pairs_are_disjoint_and_competitors_are_recorded():
    result = nominate_local_contexts(*fixture())
    assert [(row["family_a"], row["family_b"]) for row in result["contexts"]] == [(0, 1), (2, 3)]
    ranking = result["contexts"][0]["family_a_ranking"]
    assert ranking["runner_up_neighbor"] == 2 and ranking["runner_up_score"] == 4
    assert ranking["score_margin"] == 5
    assert result["counts"]["selected_contexts"] == 2
    assert not result["orthology_inferred"] and not result["families_merged"]
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("chunk", [1, 2, 3, 4, 5, 500000])
def test_chunking_duplicates_and_hit_order_do_not_change_nomination(chunk):
    arguments = fixture()
    expected = semantic(nominate_local_contexts(*arguments, chunk_size=1))
    rng = np.random.default_rng(18102026)
    order = rng.permutation(len(arguments[-1]))
    shuffled = (*arguments[:3], *[array[order] for array in arguments[3:]])
    assert semantic(nominate_local_contexts(*shuffled, chunk_size=chunk)) == expected


def test_duplicate_hits_do_not_manufacture_distinct_gene_support():
    arguments = fixture(links=())
    result = nominate_local_contexts(*arguments[:3], np.array([0, 0, 4, 4]),
        np.array([4, 4, 0, 0]), np.array([9., 9., 9., 9.]), chunk_size=1)
    assert result["counts"]["reciprocally_supported_pairs"] == 0
    assert not result["contexts"]
    assert result["pair_audit"][0]["status"] == "insufficient_reciprocal_witnesses"


def test_reciprocal_family_support_does_not_require_identical_reversed_gene_pairs():
    arguments = fixture(sizes=(3, 3), links=())
    # Each direction has two distinct queries/targets, but no exact reversed pair.
    result = nominate_local_contexts(*arguments[:3], np.array([0, 1, 3, 4]),
        np.array([4, 5, 1, 2]), np.array([9., 9., 9., 9.]))
    assert len(result["contexts"]) == 1


def test_exact_neighbor_ties_abstain_without_a_hidden_id_tiebreak():
    result = nominate_local_contexts(*fixture(sizes=(3, 3, 3), links=((0, 1, 9.), (0, 2, 9.))))
    assert not result["contexts"]
    assert result["family_ranking"][0]["best_neighbor"] is None
    assert result["family_ranking"][0]["best_tied_neighbors"] == [1, 2]


def test_joint_budget_exclusion_does_not_select_a_weaker_competitor():
    arguments = fixture(sizes=(4, 3, 4), links=((0, 2, 9.), (0, 1, 4.)))
    result = nominate_local_contexts(*arguments, max_joint_genes=7)
    assert not result["contexts"]
    assert result["family_ranking"][0]["best_neighbor"] == 2
    assert next(row for row in result["pair_audit"] if row["family_b"] == 2)["status"] == "joint_gene_budget"


def test_budget_ranking_uses_gene_content_and_survives_family_row_permutation():
    arguments = fixture(links=((0, 1, 9.), (2, 3, 9.)))
    first = nominate_local_contexts(*arguments, max_contexts=1)
    groups = arguments[0]
    permuted = ([groups[3], groups[2], groups[1], groups[0]], *arguments[1:])
    second = nominate_local_contexts(*permuted, max_contexts=1)
    genes = lambda result: sorted(result["contexts"][0]["genes_a"] + result["contexts"][0]["genes_b"])
    assert genes(first) == genes(second) == arguments[1][:6]
    assert sum(row["status"] == "context_budget" for row in first["pair_audit"]) == 1


def test_directional_maxima_are_combined_geometrically():
    arguments = fixture(sizes=(3, 3), links=((0, 1, 9.),))
    arguments[-1][:2] = 16
    result = nominate_local_contexts(*arguments)
    assert result["contexts"][0]["score"] == 12
    assert result["contexts"][0]["a_to_b_score"] == 16
    assert result["contexts"][0]["b_to_a_score"] == 9


def test_nonpositive_self_same_species_and_within_family_hits_are_ignored():
    arguments = fixture(sizes=(3, 3), links=())
    result = nominate_local_contexts(*arguments[:3], np.array([0, 0, 0, 0]),
        np.array([0, 1, 3, 4]), np.array([9., 9., 9., 0.]))
    assert result["counts"]["input_hits"] == 4
    assert result["counts"]["eligible_positive_cross_hits"] == 0
    assert not result["contexts"]


def test_readonly_hit_and_species_arrays_and_source_partition_are_unchanged():
    arguments = fixture()
    original_groups, original_names = deepcopy(arguments[0]), list(arguments[1])
    before = [array.tobytes() for array in arguments[2:]]
    for array in arguments[2:]:
        array.setflags(write=False)
    assert nominate_local_contexts(*arguments)["contexts"]
    assert arguments[0] == original_groups and arguments[1] == original_names
    assert [array.tobytes() for array in arguments[2:]] == before


def test_unsupported_families_and_empty_cases_are_explicit():
    arguments = fixture(sizes=(3, 3), links=())
    arguments[2][:3] = 0
    result = nominate_local_contexts(*arguments)
    assert result["excluded_families"] == [{"family_id": 0, "genes": 3, "species": 1, "reasons": ["too_few_species"]}]
    empty = nominate_local_contexts([], [], [], [], [], [])
    assert empty["counts"]["genes"] == 0 and empty["contexts"] == []


@pytest.mark.parametrize("key,value", [("max_contexts", 0), ("max_contexts", True),
    ("max_joint_genes", 5), ("chunk_size", 0), ("chunk_size", 1.5)])
def test_invalid_resource_bounds_rejected(key, value):
    with pytest.raises(ValueError, match="positive integer"):
        nominate_local_contexts(*fixture(), **{key: value})


@pytest.mark.parametrize("kind", ["overlap", "missing", "duplicate", "empty", "float_index", "species", "names"])
def test_invalid_complete_partition_or_ownership_rejected(kind):
    args = list(fixture())
    if kind == "overlap":
        args[0][1][0] = 0
    elif kind == "missing":
        args[0].pop()
    elif kind == "duplicate":
        args[0][0][1] = 0
    elif kind == "empty":
        args[0][0] = []
    elif kind == "float_index":
        args[0][0] = [0., 1., 2.]
    elif kind == "species":
        args[2] = args[2][:-1]
    else:
        args[1][0] = args[1][1]
    with pytest.raises(ValueError):
        nominate_local_contexts(*args)


@pytest.mark.parametrize("kind", ["negative", "too_large", "float_index", "shape", "nan", "infinity"])
def test_invalid_hit_arrays_rejected(kind):
    args = list(fixture())
    if kind == "negative":
        args[3][0] = -1
    elif kind == "too_large":
        args[4][0] = len(args[1])
    elif kind == "float_index":
        args[3] = args[3].astype(float)
    elif kind == "shape":
        args[3] = args[3][:, None]
    elif kind == "nan":
        args[5][0] = np.nan
    else:
        args[5][0] = np.inf
    with pytest.raises(ValueError):
        nominate_local_contexts(*args)


def joint_fixture(newick, ownership=None):
    pytest.importorskip("dendropy")
    from orthohmm.phylogeny import parse_gene_tree, reconcile_gene_tree
    from orthohmm.phylogeny_pipeline import FamilyOutcome, apply_satellite_membership_constraints
    left, right = ["a1", "a2", "a3"], ["b1", "b2", "b3"]
    species = (ownership or {gene: f"sp{index + 1}" for group in (left, right) for index, gene in enumerate(group)})
    tree = "((sp1,sp2),(sp3,sp4));" if "sp4" in species.values() else "((sp1,sp2),sp3);"
    result = reconcile_gene_tree(parse_gene_tree(newick), parse_gene_tree(tree), species,
        family_id="LocalFixture", root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")
    raw = FamilyOutcome("LocalFixture", tuple(sorted(left + right)), result.root_groups,
        result.ortholog_pairs, result.ortholog_pair_confidence, result, False)
    constrained, _ = apply_satellite_membership_constraints([raw], [])
    return constrained[0], left, right, species


def baseline_pairs(left, right, species):
    return [tuple(sorted((a, b))) for group in (left, right)
            for index, a in enumerate(group) for b in group[index + 1:] if species[a] != species[b]]


def test_real_reconciliation_admits_coorthology_and_preserves_existing_high_pairs():
    outcome, left, right, species = joint_fixture("((a1,b1),(a2,b2),(a3,b3));")
    original = baseline_pairs(left, right, species)
    decision = assess_joint_outcome(outcome, left, right, original, species)
    assert decision["accepted"] and not decision["reasons"]
    assert decision["baseline_high_pairs"] == decision["baseline_high_pairs_retained"] == 6
    assert len(decision["novel_cross_high_pairs"]) == 6
    assert decision["supported_species_a"] == decision["supported_species_b"] == ["sp1", "sp2", "sp3"]
    assert outcome.root_groups == (tuple(sorted(left + right)),)
    json.dumps(decision, allow_nan=False)


def test_real_ancient_duplication_does_not_become_a_rescue():
    outcome, left, right, species = joint_fixture("((a1,a2,a3),(b1,b2,b3));")
    decision = assess_joint_outcome(outcome, left, right, baseline_pairs(left, right, species), species)
    assert not decision["accepted"] and not decision["novel_cross_high_pairs"]
    assert not decision["baseline_high_pairs_lost"]
    assert decision["duplications"] >= 1


def test_novel_cross_pairs_cannot_override_existing_pair_loss():
    outcome, left, right, species = joint_fixture("((a1,b2),(a2,b1),(a3,b3));")
    decision = assess_joint_outcome(outcome, left, right, baseline_pairs(left, right, species), species)
    assert not decision["accepted"] and "baseline_high_pair_loss" in decision["reasons"]
    assert decision["baseline_high_pairs_lost"] and decision["novel_cross_high_pairs"]


def test_uncertain_mapping_is_not_new_high_confidence_orthology():
    ownership = {"a1": "sp1", "a2": "sp3", "a3": "sp1", "b1": "sp2", "b2": "sp4", "b3": "sp2"}
    outcome, left, right, species = joint_fixture("(((a1,a3),a2),((b1,b3),b2));", ownership)
    decision = assess_joint_outcome(outcome, left, right, baseline_pairs(left, right, species), species)
    assert decision["uncertain_events"] > 0
    assert not decision["accepted"] and not decision["novel_cross_high_pairs"]


def test_joint_admission_uses_post_membership_constraint_outcome():
    from orthohmm.phylogeny_pipeline import apply_satellite_membership_constraints
    outcome, left, right, species = joint_fixture("((a1,b1),(a2,b2),(a3,b3));")
    # A counterfactual changed constraint detaches a1: the gate must assess
    # the resulting loss, not the original unconstrained reconciliation.
    changed, audit = apply_satellite_membership_constraints([outcome], [
        {"source_genes": ["a1"], "target_genes": ["b1"]}])
    assert audit["detached_constraints"] == 1
    decision = assess_joint_outcome(changed[0], left, right, baseline_pairs(left, right, species), species)
    assert not decision["accepted"] and "baseline_high_pair_loss" in decision["reasons"]


@pytest.mark.parametrize("kind", ["not_reconciled", "wrong_genes", "duplicate_roots", "missing_confidence", "duplicate_baseline", "cross_baseline"])
def test_joint_admission_rejects_inconsistent_context_products(kind):
    outcome, left, right, species = joint_fixture("((a1,b1),(a2,b2),(a3,b3));")
    original = baseline_pairs(left, right, species)
    if kind == "not_reconciled":
        outcome = replace(outcome, reconciliation=None)
    elif kind == "wrong_genes":
        outcome = replace(outcome, genes=outcome.genes[:-1])
    elif kind == "duplicate_roots":
        outcome = replace(outcome, root_groups=(*outcome.root_groups, ("a1",)))
    elif kind == "missing_confidence":
        outcome = replace(outcome, ortholog_pair_confidence=outcome.ortholog_pair_confidence[:-1])
    elif kind == "duplicate_baseline":
        original.append(original[0])
    else:
        original.append(("a1", "b2"))
    with pytest.raises(ValueError):
        assess_joint_outcome(outcome, left, right, original, species)

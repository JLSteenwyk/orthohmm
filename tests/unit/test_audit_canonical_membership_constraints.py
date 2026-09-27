from itertools import permutations

import pytest

from benchmark_tools.audit_canonical_membership_constraints import semantic_constraints,compare_constraints
from orthohmm.phylogeny_pipeline import FamilyOutcome,apply_satellite_membership_constraints


def test_ignore_unused_metadata_and_gene_order_not_direction():
    a=[dict(source_genes=["a","b"],target_genes=["c"],support=1,iteration=0)]
    b=[dict(source_genes=["b","a"],target_genes=["c"],support=2,iteration=1)]
    assert compare_constraints(a,b)["semantic_sequence_equal"]
    assert not compare_constraints(a,[dict(source_genes=["c"],target_genes=["a","b"])])["semantic_multiset_equal"]


def test_counter_preserves_multiplicity_and_detects_order_changes():
    a=dict(source_genes=["a"],target_genes=["c"])
    b=dict(source_genes=["b"],target_genes=["c"])
    assert not compare_constraints([a,a],[a])["semantic_multiset_equal"]
    result=compare_constraints([a,b],[b,a])
    assert result["semantic_multiset_equal"] and not result["semantic_sequence_equal"]
    assert result["positions_changed"]==2
    assert result["historical_semantic_sha256"]==result["current_semantic_sha256"]


@pytest.mark.parametrize("source,target", [([], ["b"]),(["a","a"],["b"]),(["a"],["a"]),
    ([None],["b"]),(["a"],[""]),(["a"],[])])
def test_invalid_gene_sets_rejected(source,target):
    with pytest.raises(ValueError):
        semantic_constraints([dict(source_genes=source,target_genes=target)])


def test_actual_consumer_preserves_memberships_not_group_order():
    outcome=FamilyOutcome(family_id="F0",genes=tuple("abcde"),root_groups=(tuple("abcde"),),
        ortholog_pairs=(("a","b"),("a","c"),("a","d"),("a","e")),
        ortholog_pair_confidence=(("a","b","medium"),("a","c","high"),("a","d","low"),("a","e","high")),
        reconciliation=None,checkpoint_hit=False)
    constraints=[dict(source_genes=[g],target_genes=["a"]) for g in "bcd"]
    reference,audit=apply_satellite_membership_constraints([outcome],constraints)
    semantic={frozenset(g) for g in reference[0].root_groups}
    ordered=set()
    for rows in permutations(constraints):
        changed=[dict(r,support=-1e100,margin=0,iteration=99,source_cluster=999) for r in rows]
        actual,summary=apply_satellite_membership_constraints([outcome],changed)
        assert {frozenset(g) for g in actual[0].root_groups}==semantic
        assert actual[0].ortholog_pairs==reference[0].ortholog_pairs
        assert actual[0].ortholog_pair_confidence==reference[0].ortholog_pair_confidence
        assert summary==audit
        ordered.add(actual[0].root_groups)
    assert len(ordered)==2

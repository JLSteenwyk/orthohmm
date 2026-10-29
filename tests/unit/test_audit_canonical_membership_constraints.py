from itertools import permutations

import pytest

from orthohmm.phylogeny_pipeline import FamilyOutcome,apply_satellite_membership_constraints








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

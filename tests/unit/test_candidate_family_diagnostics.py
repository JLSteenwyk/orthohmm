import pytest

from benchmark_tools.candidate_family_diagnostics import (
    audit_candidate_partition,
)


def test_candidate_audit_reports_recall_growth_and_contamination():
    seeds = [{"a", "x"}, {"b"}, {"c"}, {"d"}]
    candidates = [{"a", "x", "b"}, {"c", "d"}]
    refogs = [{"a", "b", "c"}, {"d"}]

    result = audit_candidate_partition(seeds, candidates, refogs)

    assert result["seed_merges"] == 2
    assert result["candidates_combining_seeds"] == 2
    assert result["candidate_family_sizes"]["maximum"] == 3
    assert result["reference"]["possible_pairs"] == 3
    assert result["reference"]["recovered_pairs"] == 1
    assert result["reference"]["pair_recall_micro"] == pytest.approx(1 / 3)
    assert result["reference"]["cross_refog_pairs"] == 1
    assert result["reference"]["pair_precision_among_reference_genes"] == 0.5
    assert result["reference"]["reference_to_nonreference_relations"] == 2
    assert result["refog_records"][0]["seed_components"] == 3
    assert result["refog_records"][0]["candidate_components"] == 2


def test_candidate_audit_rejects_seed_splitting():
    with pytest.raises(ValueError, match="splits immutable seed family"):
        audit_candidate_partition(
            [{"a", "b"}],
            [{"a"}, {"b"}],
            [{"a", "b"}],
        )


def test_candidate_audit_rejects_gene_universe_changes():
    with pytest.raises(ValueError, match="gene universes differ"):
        audit_candidate_partition(
            [{"a"}, {"b"}],
            [{"a"}, {"c"}],
            [{"a", "b"}],
        )

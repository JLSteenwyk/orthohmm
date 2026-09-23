import pytest

from benchmark_tools.inventory_swisstree_entry_aliases import aliases


def annotation(gene, name):
    return dict(accession=gene, selection=dict(name=name), selection_class="baseline_release",
                sequence_sha256="a", taxid="9606")


def test_exact_name_adds_candidate_without_normalization():
    result = aliases(["A_HUMAN", "B_CANFA"], {"a", "b"},
                     {"a": annotation("a", "A_HUMAN"), "b": annotation("b", "B_CANLF")})
    assert result["unique_candidates"] == 1
    assert result["missing"] == 1
    assert result["genes"]["a"]["candidates"] == {"A_HUMAN": ["historical_entry_name_exact"]}
    assert result["mapping_admitted"] is False


def test_multiple_matching_leaves_remain_ambiguous():
    result = aliases(["a", "A_HUMAN"], {"a"}, {"a": annotation("a", "A_HUMAN")})
    assert result["unique_candidates"] == 0
    assert result["genes"]["a"]["status"] == "ambiguous"


def test_shared_leaf_candidate_not_unique():
    result = aliases(["SAME"], {"a", "b"}, {g: annotation(g, "SAME") for g in ("a", "b")})
    assert result["unique_candidates"] == 0
    assert result["shared_leaf_candidates"] == {"SAME": ["a", "b"]}


def test_changed_accession_rejected():
    with pytest.raises(ValueError):
        aliases(["A"], {"a"}, {"a": annotation("b", "A")})

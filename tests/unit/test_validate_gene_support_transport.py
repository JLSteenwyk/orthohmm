import pytest

from benchmark_tools.validate_gene_support_transport import (
    independent_tree_inventory, matching_lengths, validate_family,
)
from orthohmm.phylogeny import parse_gene_tree
from orthohmm.phylogeny_pipeline import _root_external_gene_tree


def test_independent_inventory_reads_percentages_and_missing_support():
    inventory = independent_tree_inventory(parse_gene_tree("((a:1,b:1)31:1,c:1,(d:1,e:1)nan:1);"))
    assert inventory[0] == set("abcde")
    assert list(inventory[2].values()) == [{.31}]
    assert len(inventory[1]) == 7


@pytest.mark.parametrize("raw", ["(a,b,c);", "(a:1,b:inf,c:1);"])
def test_inventory_rejects_unavailable_nonfinite_lengths(raw):
    with pytest.raises(ValueError, match="branch length"):
        independent_tree_inventory(parse_gene_tree(raw))


def test_matching_lengths_requires_exact_topology_and_tolerates_serialization_roundoff():
    assert matching_lengths({("a",): 1.}, {("a",): 1. + 1e-12})
    assert matching_lengths({("a",): None}, {("a",): None})
    assert not matching_lengths({("a",): 1.}, {("b",): 1.})
    assert not matching_lengths({("a",): 1.}, {("a",): 1.01})
    assert not matching_lengths({("a",): None}, {("a",): 0.})


def fixture():
    raw = "((g00000000:0.1,g00000001:0.2)0.31:0.5,g00000002:0.2,(g00000003:0.1,g00000004:0.1)0.72:0.4);"
    species = parse_gene_tree("(a:1,(b:1,(c:1,(d:1,e:1):1):1):1);")
    owners = {gene:gene for gene in "abcde"}
    tokens = {f"g{index:08d}": gene for index, gene in enumerate("abcde")}
    saved = _root_external_gene_tree(raw, tokens, species, owners).as_string(schema="newick")
    return raw, saved, species, owners


def test_validate_family_repeats_native_rooting_and_preserves_support_and_current_outputs():
    raw, saved, species, owners = fixture()
    result = validate_family(raw, saved, tuple(owners), species, owners)
    assert result["legacy_support_map_changed"]
    assert result["repaired_support_map_preserved"]
    assert result["rooted_topology_and_lengths_unchanged"]
    assert result["current_pairs_confidence_and_groups_unchanged"]


def test_validate_family_rejects_changed_saved_support():
    raw, saved, species, owners = fixture()
    saved = saved.replace("0.72", "0.73")
    with pytest.raises(ValueError, match="Legacy replay differs"):
        validate_family(raw, saved, tuple(owners), species, owners)


def test_validate_family_rejects_changed_saved_lengths():
    raw, saved, species, owners = fixture()
    saved = saved.replace(":0.4", ":0.41")
    with pytest.raises(ValueError, match="Legacy replay differs"):
        validate_family(raw, saved, tuple(owners), species, owners)

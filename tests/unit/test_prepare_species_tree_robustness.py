import io

from Bio import Phylo
import pytest

from benchmark_tools.prepare_species_tree_robustness import digest, neighbors, panel, serialize, topology


def tree(text="(((a:1,b:2):3,(c:4,d:5):6):7,((e:8,f:9):10,(g:11,h:12):13):14):0;"):
    return Phylo.read(io.StringIO(text), "newick")


def test_neighbors_change_one_rooted_clade_and_preserve_original():
    source = tree()
    before = serialize(source)
    original = topology(source)
    variants = neighbors(source)
    assert len(variants) == 2 * (8 - 2)
    assert all(len(key ^ original) == 2 for key in variants)
    assert all(sorted(n.name for n in v.get_terminals()) == list("abcdefgh") for v in variants.values())
    assert serialize(source) == before


def test_seven_fixed_variants_with_distinct_distances_and_roundtrip():
    source = tree()
    variants = panel(source)
    assert [row[2] for row in variants] == [0, 2, 2, 2, 4, 4, 4]
    keys = [topology(row[1]) for row in variants]
    assert len(set(keys)) == 7
    assert keys == [topology(row[1]) for row in panel(source)]
    lengths = sorted(n.branch_length for n in source.find_clades())
    for _, candidate, distance, _ in variants:
        text = serialize(candidate)
        assert text.startswith("[&R] ")
        restored = Phylo.read(io.StringIO(text), "newick")
        assert topology(restored) == topology(candidate)
        assert len(topology(restored) ^ topology(source)) == distance
        assert sorted(n.branch_length for n in restored.find_clades()) == lengths


def test_hash_ignores_sibling_order():
    assert digest(topology(tree("((a,b),(c,d));"))) == digest(topology(tree("((d,c),(b,a));")))


@pytest.mark.parametrize("text", ["((a,a),(c,d));", "(a,b,c,d);", "((a:-1,b),(c,d));", "(a,(b,c));"])
def test_invalid_input_rejected(text):
    with pytest.raises(ValueError):
        topology(tree(text))

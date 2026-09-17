from io import StringIO

from Bio import Phylo
import pytest

from benchmark_tools.prepare_simulation_tree_controls import restrict_tree, controls, topology, prepare


def tree():
    return Phylo.read(StringIO("(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);"), "newick")


def test_pruned_generating_tree_and_deterministic_controls():
    original = tree()
    before = topology(original)
    subset = restrict_tree(original, ["a", "b", "c", "e", "g", "h"])
    assert topology(original) == before
    rows = controls(subset)
    assert [row[0] for row in rows] == ["generating", "nni1", "nni2"]
    assert [row[2] for row in rows] == [0, 2, 4]
    for _, candidate, distance, _ in rows:
        assert {tip.name for tip in candidate.get_terminals()} == {"a", "b", "c", "e", "g", "h"}
        assert len(topology(candidate) ^ topology(subset)) == distance
    assert [topology(r[1]) for r in rows] == [topology(r[1]) for r in controls(subset)]


@pytest.mark.parametrize("names", [["a", "b", "c"], ["a", "b", "c", "x"], ["a", "b", "c", "c"]])
def test_invalid_retained_taxa(names):
    with pytest.raises(ValueError):
        restrict_tree(tree(), names)


def test_root_collapse_preserves_retained_distances():
    original = tree()
    subset = restrict_tree(original, ["a", "b", "c", "d"])
    assert subset.distance("a", "c") == original.distance("a", "c")
    assert len(subset.get_terminals()) == 4


def test_existing_output_not_reused(tmp_path):
    with pytest.raises(FileExistsError):
        prepare(tmp_path, tmp_path)

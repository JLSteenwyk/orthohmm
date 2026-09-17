import io

from Bio import Phylo
import pytest

from benchmark_tools.prepare_portable_simulation_trees import portable_text, branch_map, prepare


def test_plain_newick_preserves_root_topology_lengths_and_original():
    tree = Phylo.read(io.StringIO("[&R] ((a:0.1,b:0.2):0.3,(c:0.4,d:0.5):0.6):0.7;"), "newick")
    before = branch_map(tree)
    text = portable_text(tree)
    assert text.startswith("(") and "[&R]" not in text
    assert branch_map(Phylo.read(io.StringIO(text), "newick")) == before
    assert branch_map(tree) == before


def test_polytomy_rejected():
    with pytest.raises(ValueError, match="bifurcating"):
        portable_text(Phylo.read(io.StringIO("(a,b,c,d);"), "newick"))


def test_existing_output_rejected(tmp_path):
    with pytest.raises(FileExistsError):
        prepare(tmp_path, tmp_path)

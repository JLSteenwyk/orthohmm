from io import StringIO

import pytest

from benchmark_tools.proteinortho_to_pairwise import iter_pairs, write_pairs


def _header() -> str:
    return (
        "# file_a\tfile_b\n"
        "# a\tb\tevalue_ab\tbitscore_ab\tevalue_ba\tbitscore_ba\n"
        "# species_a.fasta\tspecies_b.fasta\n"
        "# 1e-10\t50\t1e-10\t50\n"
    )


def test_converts_proteinortho_graph(tmp_path):
    graph = tmp_path / "proteinortho-graph"
    graph.write_text(
        _header()
        + "sp|B|GENE_B\ttr|A|GENE_A\t1e-20\t80\t2e-20\t75\n"
    )
    output = StringIO()

    assert write_pairs(graph, output) == 1
    assert output.getvalue() == "A\tB\n"


def test_rejects_duplicate_pair_within_species_section(tmp_path):
    graph = tmp_path / "proteinortho-graph"
    relation = "sp|B|GENE_B\ttr|A|GENE_A\t1e-20\t80\t2e-20\t75\n"
    graph.write_text(_header() + relation + relation)

    with pytest.raises(ValueError, match="Duplicate ProteinOrtho pair"):
        list(iter_pairs(graph))


def test_rejects_duplicate_species_section(tmp_path):
    graph = tmp_path / "proteinortho-graph"
    graph.write_text(
        "# species_a.fasta\tspecies_b.fasta\n"
        "a\tb\t1e-10\t50\t1e-10\t50\n"
        "# species_b.fasta\tspecies_a.fasta\n"
        "b\ta\t1e-10\t50\t1e-10\t50\n"
    )

    with pytest.raises(ValueError, match="Duplicate ProteinOrtho species section"):
        list(iter_pairs(graph))


def test_rejects_incomplete_species_pair_matrix(tmp_path):
    graph = tmp_path / "proteinortho-graph"
    graph.write_text(
        "# species_a.fasta\tspecies_b.fasta\n"
        "a\tb\t1e-10\t50\t1e-10\t50\n"
        "# species_a.fasta\tspecies_c.fasta\n"
        "a\tc\t1e-10\t50\t1e-10\t50\n"
    )

    with pytest.raises(ValueError, match="Incomplete ProteinOrtho"):
        list(iter_pairs(graph))

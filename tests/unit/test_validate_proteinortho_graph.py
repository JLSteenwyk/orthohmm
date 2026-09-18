import pytest

from benchmark_tools.validate_proteinortho_graph import validate_graph


def setup_graph(tmp_path, text):
    inputs = []
    for species, accession in (("s1", "a"), ("s2", "b"), ("s3", "c")):
        path = tmp_path / f"{species}.fasta"
        path.write_text(f">sp|{accession}|name\nACDE\n")
        inputs.append(path)
    graph = tmp_path / "qfo.proteinortho-graph"
    graph.write_text(text)
    return graph, inputs


VALID = ("# file_a\tfile_b\n# s1.fasta\ts2.fasta\n"
         "sp|a|name\tsp|b|name\t0\t50\t1e-20\t60\n"
         "# s1.fasta\ts3.fasta\n# s2.fasta\ts3.fasta\n")


def test_valid_including_empty_sections(tmp_path):
    assert validate_graph(*setup_graph(tmp_path, VALID)) == {
        "input_species": 3, "input_accessions": 3,
        "species_sections": 3, "validated_pair_rows": 1}


@pytest.mark.parametrize("text", [
    VALID.replace("sp|a|name", "sp|c|name"),
    VALID.replace("sp|a|name", "sp|foreign|name"),
    VALID.replace("\t50\t", "\tnan\t"),
    VALID.replace("\t50\t", "\tinf\t"),
    VALID.replace("# s1.fasta\ts3.fasta\n# s2.fasta\ts3.fasta\n", ""),
    VALID.replace("s3.fasta", "unknown.fasta"),
    VALID + "# s2.fasta\ts1.fasta\n",
    VALID.replace("sp|a|name", "sp|a|bad name"),
    "", "a\tb\t0\t50\t0\t50\n",
])
def test_reject_invalid(tmp_path, text):
    with pytest.raises(ValueError):
        validate_graph(*setup_graph(tmp_path, text))


def test_reject_search_graph(tmp_path):
    graph, inputs = setup_graph(tmp_path, VALID)
    target = graph.with_name("qfo.blast-graph")
    graph.rename(target)
    with pytest.raises(ValueError, match="post-clustering"):
        validate_graph(target, inputs)


def test_existing_duplicate_check_retained(tmp_path):
    text = VALID.replace("# s1.fasta\ts3.fasta", "sp|a|name\tsp|b|name\t0\t50\t0\t50\n# s1.fasta\ts3.fasta")
    with pytest.raises(ValueError, match="Duplicate ProteinOrtho pair"):
        validate_graph(*setup_graph(tmp_path, text))

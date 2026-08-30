import pytest

from benchmark_tools.native_pairs_to_qfo import convert, strip_to_uniprot


def test_strip_to_uniprot_handles_qfo_headers():
    assert strip_to_uniprot("tr|A0A024R1R8|A0A024R1R8_HUMAN") == "A0A024R1R8"
    assert strip_to_uniprot("plain") == "plain"


def test_convert_native_pairs_selects_gene_columns_and_canonicalizes(tmp_path):
    source = tmp_path / "native.tsv"
    source.write_text(
        "gene_a\tspecies_a\tgene_b\tspecies_b\n"
        "tr|Z9|Z_SP\tsp1\tsp|A1|A_SP\tsp2\n"
    )
    output = tmp_path / "qfo.tsv"

    assert convert(source, output) == 1
    assert output.read_text() == "A1\tZ9\n"


def test_convert_native_pairs_rejects_wrong_schema(tmp_path):
    source = tmp_path / "native.tsv"
    source.write_text("left\tright\n")

    with pytest.raises(ValueError, match="unexpected native pair header"):
        convert(source, tmp_path / "qfo.tsv")

import pytest

from benchmark_tools.audit_orthofinder_input_parity import read_maps, compare_species


def test_maps_and_exact_or_changed_sequences(tmp_path):
    species, identifiers = tmp_path / "species", tmp_path / "ids"
    species.write_text("0: input.fasta\n")
    identifiers.write_text("0_0: sp|A|ONE protein: description\n0_1: tr|B|TWO other\n")
    names, mapping = read_maps(species, identifiers)
    assert names == {"0": "input.fasta"}
    original, native = tmp_path / "input.fasta", tmp_path / "Species0.fa"
    original.write_text(">sp|A|ONE\nACD\n>tr|B|TWO\nAAA\n")
    native.write_text(">0_0\nACD\n>0_1\nAAX\n")
    result = compare_species(original, native, mapping["0"])
    assert result["sequences"] == 2 and result["identical_sequences"] == 1
    assert result["differences"][0]["original_id"] == "tr|B|TWO"
    native.write_text(">0_0\nACD\n")
    with pytest.raises(ValueError, match="Incomplete"):
        compare_species(original, native, mapping["0"])


@pytest.mark.parametrize("species_text,ids_text", [
    ("0: a\n0: b\n", "0_0: A\n"), ("0: ../a\n", "0_0: A\n"),
    ("0: a\n", "1_0: A\n"), ("0: a\n", "0_0: A\n0_0: B\n"),
    ("0: a\n", "0_0: A\n0_1: A\n"), ("0: a\n", "bad: A\n"),
])
def test_invalid_maps(tmp_path, species_text, ids_text):
    species, identifiers = tmp_path / "species", tmp_path / "ids"
    species.write_text(species_text)
    identifiers.write_text(ids_text)
    with pytest.raises(ValueError):
        read_maps(species, identifiers)


@pytest.mark.parametrize("text", [">0_0\nACD\n>0_0\nACD\n", ">0_1\nACD\n"])
def test_bad_internal_ids(tmp_path, text):
    original, native = tmp_path / "original", tmp_path / "native"
    original.write_text(">A\nACD\n")
    native.write_text(text)
    with pytest.raises(ValueError):
        compare_species(original, native, {"0_0": "A"})

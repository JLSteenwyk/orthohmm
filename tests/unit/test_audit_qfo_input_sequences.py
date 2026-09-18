import pytest

from benchmark_tools.audit_qfo_input_sequences import compare_database, parse_entry, sequence_identity


def test_field_parsing():
    species, aliases, sequence = parse_entry("<E><OS>A</OS><MAPIDS>X; Y</MAPIDS><SEQ>ACD</SEQ></E>")
    assert species == "A" and aliases == {"X", "Y"}
    assert sequence == sequence_identity("ACD")
    assert sequence != sequence_identity("acd")


@pytest.mark.parametrize("line", ["<E/>", "<X/>", "<E><OS>A</OS><MAPIDS>X; X</MAPIDS><SEQ>A</SEQ></E>",
    "<E><OS>A</OS><MAPIDS>X</MAPIDS><SEQ>A A</SEQ></E>",
    "<E><OS>A</OS><MAPIDS>X</MAPIDS><SEQ>A</SEQ><SEQ>B</SEQ></E>"])
def test_invalid_entry(line):
    with pytest.raises(ValueError):
        parse_entry(line)


def test_same_changed_missing_and_coverage(tmp_path):
    path = tmp_path / "db"
    path.write_text("".join(f"<E><OS>A</OS><MAPIDS>{a}</MAPIDS><SEQ>{s}</SEQ></E>\n" for a, s in [("X", "ACD"), ("Y", "AAA"), ("Z", "CCC")]))
    inputs = {1: {"accession": "X", "sequence": sequence_identity("ACD")}, 2: {"accession": "Y", "sequence": sequence_identity("AAT")}}
    result = compare_database(path, inputs, 3)
    assert result["sequence_identical"] == result["sequence_different"] == 1
    assert result["native_entries_without_input_by_species"] == {"A": 1}
    assert result["differences"][0]["numeric_id"] == 2
    with pytest.raises(ValueError, match="count"):
        compare_database(path, inputs, 4)
    inputs[1]["accession"] = "WRONG"
    with pytest.raises(ValueError, match="alias"):
        compare_database(path, inputs, 3)

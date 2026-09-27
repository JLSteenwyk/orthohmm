import pytest

from benchmark_tools.audit_ob_orthofinder_provenance import canonical, compare_species, fasta, log_command, read_groups, species_ids


def test_sequence_content_and_exact_id_mapping():
    mapping = {"0_0": "a", "0_1": "b", "1_0": "c"}
    result = compare_species(dict(a="AC", b="DF"), {"0_0": "AC", "0_1": "DF"}, mapping, 0)
    assert result["mapping_ids_match_processed"] and result["original_ids_match"]
    assert result["sequence_mismatch_ids"] == []
    assert compare_species(dict(a="AC", b="DF"), {"0_0": "AA", "0_1": "DF"}, mapping, 0)["sequence_mismatch_ids"] == ["a"]
    assert not compare_species(dict(a="AC", b="DF"), {"0_0": "AC"}, mapping, 0)["original_ids_match"]
    assert not compare_species(dict(a="AC", b="DF"), {"0_0": "AC", "0_1": "DF", "extra": "AC"}, mapping, 0)["mapping_ids_match_processed"]


def test_duplicate_mapped_names_rejected():
    with pytest.raises(ValueError, match="Duplicate"):
        compare_species(dict(a="AC"), {"0_0": "AC", "0_1": "AC"}, {"0_0": "a", "0_1": "a"}, 0)


def test_fasta_duplicates_and_empty(tmp_path):
    path = tmp_path / "input.fa"
    path.write_text(">a description\nAC\nD\n")
    assert fasta(path) == dict(a="ACD")
    path.write_text(">a\nAC\n>a\nDF\n")
    with pytest.raises(ValueError, match="Duplicate"):
        fasta(path)
    path.write_text("")
    with pytest.raises(ValueError, match="Empty"):
        fasta(path)


@pytest.mark.parametrize("text", ["0: a.fa\n0: b.fa\n", "1: a.fa\n", "0: ../a.fa\n", "0: a.fa\n1: a.fa\n", "bad\n", ""])
def test_invalid_species_mapping(tmp_path, text):
    path = tmp_path / "SpeciesIDs.txt"
    path.write_text(text)
    with pytest.raises(ValueError):
        species_ids(path)


def test_valid_species_mapping(tmp_path):
    path = tmp_path / "SpeciesIDs.txt"
    path.write_text("0: a.fa\n1: b.fa\n")
    assert species_ids(path) == {0: "a.fa", 1: "b.fa"}


def test_partition_is_label_order_independent_but_complete():
    assert canonical([["b", "a"], ["c"]], {"a", "b", "c"}) == [("a", "b"), ("c",)]
    for groups in ([["a"], ["a", "b", "c"]], [["a", "b"]], [[], ["a", "b", "c"]]):
        with pytest.raises(ValueError, match="exact gene universe"):
            canonical(groups, {"a", "b", "c"})


def test_command_quotes_and_duplicate_detection():
    assert log_command('Command Line: tool -f "input with spaces"\n', "Command Line: ") == ["tool", "-f", "input with spaces"]
    assert log_command('\tCommand being timed: "tool -t 32"\n', "Command being timed: ") == ["tool", "-t", "32"]
    for text in ("", "Command Line: a\nCommand Line: b\n"):
        with pytest.raises(ValueError, match="one command"):
            log_command(text, "Command Line: ")


def test_group_parser_preserves_duplicates_for_validation(tmp_path):
    path = tmp_path / "groups.txt"
    path.write_text("OG0: a a\nOG1: b\n")
    with pytest.raises(ValueError, match="exact gene universe"):
        canonical(read_groups(path, True), {"a", "b"})
    path.write_text("OG0: a\nOG0: b\n")
    with pytest.raises(ValueError, match="duplicate native group label"):
        read_groups(path, True)
    path.write_text("a b\nc\n")
    assert read_groups(path, False) == [["a", "b"], ["c"]]

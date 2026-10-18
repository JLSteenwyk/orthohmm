import pytest

from benchmark_tools.native_pair_turnover import native_pair_turnover


HEADER = "gene_a\tspecies_a\tgene_b\tspecies_b\tconfidence\n"


def test_pair_addition_loss_confidence_changes_and_order_independence(tmp_path):
    old, new = tmp_path / "old.tsv", tmp_path / "new.tsv"
    old.write_text(HEADER + "a\tA\tb\tB\thigh\na\tA\tc\tC\thigh\nb\tB\td\tD\tlow\n")
    new.write_text(HEADER + "b\tB\ta\tA\tmedium\nb\tB\td\tD\thigh\na\tA\td\tD\thigh\n")
    actual = native_pair_turnover(old, new)
    assert actual["baseline_pairs"] == actual["candidate_pairs"] == 3
    assert actual["retained_pairs"] == 2 and actual["added_pairs"] == actual["lost_pairs"] == 1
    assert actual["baseline_high_pairs_absent"] == actual["baseline_high_pairs_downgraded"] == 1
    assert actual["retained_confidence_transitions"] == {"high->medium": 1, "low->high": 1}


@pytest.mark.parametrize("row", ["a\tA\ta\tA\thigh\n", "a\tA\tb\tA\thigh\n", "a\tA\tb\tB\tbad\n",
    "a\tA\tb\tB\thigh\na\tA\tb\tB\tlow\n", "a\tOther\tb\tB\thigh\n"])
def test_invalid_duplicate_or_changed_species_rows_are_not_silently_counted(tmp_path, row):
    old, new = tmp_path / "old.tsv", tmp_path / "new.tsv"
    old.write_text(HEADER + "a\tA\tb\tB\thigh\n")
    new.write_text(HEADER + row)
    with pytest.raises(ValueError): native_pair_turnover(old, new)


def test_empty_native_files_remain_explicit_zero_population(tmp_path):
    path = tmp_path / "empty.tsv"
    path.write_text(HEADER)
    result = native_pair_turnover(path, path)
    assert result["baseline_pairs"] == result["candidate_pairs"] == result["retained_pairs"] == 0

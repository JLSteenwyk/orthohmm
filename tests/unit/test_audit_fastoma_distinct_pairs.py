import pytest

from benchmark_tools.audit_fastoma_distinct_pairs import compare_sets, audit


def test_independent_sort_compares_full_contents(tmp_path):
    actual, retained, sorted_path = [tmp_path / name for name in ("actual", "retained", "sorted")]
    actual.write_text("a\tb\nc\td\n")
    retained.write_text("c\td\na\tb\na\tb\n")
    compare_sets(actual, retained, sorted_path)
    assert sorted_path.read_bytes() == actual.read_bytes()


def test_different_pair_sets_fail(tmp_path):
    actual, retained, sorted_path = [tmp_path / name for name in ("actual", "retained", "sorted")]
    actual.write_text("a\tb\n")
    retained.write_text("a\tc\n")
    with pytest.raises(ValueError, match="differs"):
        compare_sets(actual, retained, sorted_path)


def test_sort_output_not_overwritten(tmp_path):
    existing = tmp_path / "existing"
    existing.write_text("preserve")
    with pytest.raises(FileExistsError):
        compare_sets(tmp_path / "absent", tmp_path / "absent", existing)
    assert existing.read_text() == "preserve"


@pytest.mark.parametrize("existing", ["work", "report"])
def test_audit_requires_fresh_paths(tmp_path, existing):
    (tmp_path / existing).touch()
    with pytest.raises(FileExistsError):
        audit(tmp_path, tmp_path / "work", tmp_path / "report")

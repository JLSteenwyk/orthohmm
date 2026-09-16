import pytest

from benchmark_tools.audit_historical_profile_ablation import read_partition, verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance


def test_verified_partition(tmp_path):
    path = tmp_path / "groups.txt"
    path.write_text("a b\nc\n")
    expected = file_provenance(path)
    assert verify_file(path, expected) == expected
    assert read_partition(path, {"a", "b", "c"}) == [{"a", "b"}, {"c"}]
    path.write_text("a c\nb\n")
    with pytest.raises(ValueError, match="Changed"):
        verify_file(path, expected)


@pytest.mark.parametrize("text", ["a a\nb c\n", "a b\na c\n", "a b\n", "a b c unknown\n"])
def test_invalid_partition(tmp_path, text):
    path = tmp_path / "groups.txt"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_partition(path, {"a", "b", "c"})

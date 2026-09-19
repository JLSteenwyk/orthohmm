import pytest

from benchmark_tools.verify_cpu_wheel_install import partition


@pytest.mark.parametrize("text,genes", [("OG0: a b\nOG1: c\n", ["a", "b", "c"]),
                                        ("OG0: c b a\n", ["a", "b", "c"])])
def test_exact_partition(tmp_path, text, genes):
    path = tmp_path / "groups.txt"
    path.write_text(text)
    result = partition(path, genes)
    assert result["genes"] == 3
    assert result["groups"] == len(text.splitlines())


@pytest.mark.parametrize("text,genes", [("OG0: a a\n", ["a", "b"]), ("OG0: a\n", ["a", "b"]),
    ("OG0: a b c\n", ["a", "b"]), ("OG0: a a\n", ["a", "a"]), ("OG0: \n", ["a"]), (": a\n", ["a"])])
def test_invalid_partition_rejected(tmp_path, text, genes):
    path = tmp_path / "groups.txt"
    path.write_text(text)
    with pytest.raises(ValueError):
        partition(path, genes)

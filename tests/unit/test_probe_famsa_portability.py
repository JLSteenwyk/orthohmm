import pytest

from benchmark_tools.probe_famsa_portability import fixtures, validate_alignment


def test_accepts_reordered_alignment():
    assert validate_alignment(">b\nAC-\n>a\nA-C\n", {"a": "AC", "b": "AC"}) == {"b": "AC-", "a": "A-C"}


@pytest.mark.parametrize("text", [">a\nAC\n>a\nAC\n", ">b\nAC\n", ">a\nAC\n>b\nA\n",
                                  ">a\nAC\n>b\nAT\n"])
def test_rejects_invalid_alignment(text):
    with pytest.raises(ValueError):
        validate_alignment(text, {"a": "AC", "b": "AC"})


def test_fixture_scope():
    cases = fixtures()
    assert set(cases) == {"duplicates", "indels", "ambiguous"}
    assert cases["duplicates"]["a"] == cases["duplicates"]["b"]
    assert len({len(s) for s in cases["indels"].values()}) > 1
    assert "X" in cases["ambiguous"]["c"]

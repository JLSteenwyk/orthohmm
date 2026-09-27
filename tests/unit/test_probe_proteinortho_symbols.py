import pytest

from benchmark_tools.probe_proteinortho_symbols import sequences, validate


def test_fixture_sequences_differ_only_in_asterisk_location():
    values = sequences()
    assert len(values["clean"]) == 80
    for case in ("internal_asterisk", "terminal_asterisk"):
        assert values[case].replace("*", "") == values["clean"]
        assert values[case].count("*") == 1
    assert not values["internal_asterisk"].endswith("*")
    assert values["terminal_asterisk"].endswith("*")


@pytest.mark.parametrize("case,code,log,unchanged,expected", [
    ("clean", 0, "ok", True, True), ("clean", 1, "error", True, False),
    ("clean", 0, "ok", False, False),
    ("internal_asterisk", 1, "Invalid symbol; sanitze with", True, True),
    ("terminal_asterisk", 1, "Invalid symbol; sanitze with", True, True),
    ("internal_asterisk", 0, "Invalid symbol; sanitze with", True, False),
    ("internal_asterisk", 1, "unrelated error", True, False),
    ("terminal_asterisk", 1, "Invalid symbol; sanitze with", False, False),
])
def test_policy_outcome_requires_specific_failure_and_unchanged_input(case, code, log, unchanged, expected):
    assert validate(case, code, log, unchanged) is expected

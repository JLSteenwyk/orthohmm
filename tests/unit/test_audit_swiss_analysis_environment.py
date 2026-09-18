import pytest

from benchmark_tools.audit_swiss_analysis_environment import check_pins


def test_canonical_exact_pins():
    assert check_pins("# isolated\nPillow==12.3.0\n", {"pillow": "12.3.0"}) == {
        "package": [{"name": "pillow", "version": "12.3.0"}]}


@pytest.mark.parametrize("pins", ["pillow>=12.3", "pillow==12.*", "pillow==12.3.0\nPillow==12.3.0",
                                 "pillow==12.3.0; python_version>='3.10'", "pillow[extra]==12.3.0", ""])
def test_inexact_or_ambiguous_pins_rejected(pins):
    with pytest.raises(ValueError):
        check_pins(pins, {"pillow": "12.3.0"})


def test_unexpected_or_changed_packages_rejected():
    for installed in ({"pillow": "12.2.0"}, {"pillow": "12.3.0", "extra": "1"}):
        with pytest.raises(ValueError, match="inventory"):
            check_pins("pillow==12.3.0", installed)

import pytest

from tests.integration.output_checks import check_outputs


@pytest.mark.integration
@pytest.mark.parametrize("reverse", [False, True])
def test_long_names_all_outputs(native_case, reverse):
    assert check_outputs(*native_case("long_names", reverse)) == (988, 9)

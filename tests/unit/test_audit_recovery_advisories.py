import pytest

from benchmark_tools.audit_recovery_advisories import verify_inventory


def test_normalized_distribution_identity():
    installed = dict(install=[dict(metadata=dict(name="Example_Pkg", version="1.2"))])
    assert verify_inventory(installed, [dict(name="example-pkg", version="1.2")]) == dict(
        package=[dict(name="example-pkg", version="1.2")])


@pytest.mark.parametrize("actual", [[], [dict(name="x", version="2")],
    [dict(name="x", version="1"), dict(name="x", version="1")],
    [dict(name="x", version="1"), dict(name="y", version="1")]])
def test_missing_changed_duplicate_extra_inventory(actual):
    installed = dict(install=[dict(metadata=dict(name="x", version="1"))])
    with pytest.raises(ValueError):
        verify_inventory(installed, actual)

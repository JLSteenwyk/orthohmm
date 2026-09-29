import pytest

from benchmark_tools.build_timing_controller import SELECTED, build, require_inventory


def test_exact_inventory():
    require_inventory({k.lower(): v for k, v in SELECTED.items()})


@pytest.mark.parametrize("kind", ["extra", "missing", "version"])
def test_drift_rejected(kind):
    actual = dict(SELECTED)
    if kind == "extra":
        actual["unexpected"] = "1"
    elif kind == "missing":
        actual.pop("psutil")
    else:
        actual["psutil"] = "0"
    with pytest.raises(ValueError, match="inventory"):
        require_inventory(actual)


def test_existing_output_rejected(tmp_path):
    with pytest.raises(FileExistsError):
        build(tmp_path, tmp_path)

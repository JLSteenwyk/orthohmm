import pytest

from benchmark_tools.probe_verified_threadripper_fixture import METHODS, select_fixture


@pytest.mark.parametrize("index,method", enumerate(METHODS))
def test_selects_frozen_method_without_changes(index, method):
    rows = [dict(index=i, native_method=m, native_argv=[m, "frozen-flag"]) for i, m in enumerate(METHODS)]
    assert select_fixture(dict(runs=rows), method) is rows[index]


def test_wrong_method_rejected():
    with pytest.raises(ValueError):
        select_fixture(dict(runs=[]), "unknown")


def test_reordered_block_rejected():
    rows = [dict(index=i, native_method=m) for i, m in enumerate(reversed(METHODS))]
    with pytest.raises(ValueError):
        select_fixture(dict(runs=rows), METHODS[0])

import pytest

from benchmark_tools.inspect_qfo_public_forks import candidates


def test_case_insensitive_file_leads_only():
    rows = [dict(type="blob", path=p) for p in ("treefam2reference.txt", "TreeFam-A.tar.gz", "x.NHX", "other")]
    rows.append(dict(type="tree", path="treefam"))
    assert candidates(dict(truncated=False, tree=rows)) == rows[:3]


@pytest.mark.parametrize("value", [{}, dict(truncated=True, tree=[]), dict(truncated=False, tree=None)])
def test_reject_incomplete_inventory(value):
    with pytest.raises(ValueError):
        candidates(value)

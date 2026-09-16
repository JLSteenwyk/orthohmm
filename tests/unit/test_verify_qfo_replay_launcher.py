import pytest

from benchmark_tools.verify_qfo_replay_launcher import compare_trees


@pytest.fixture
def roots(tmp_path):
    roots = [tmp_path / name for name in ("frozen", "launcher")]
    for root in roots:
        (root / "orthohmm").mkdir(parents=True)
        (root / "orthohmm/core.py").write_text("x = 1\n")
        (root / "orthohmm/kernel.so").write_bytes(b"native")
    return roots


def test_exact_source_and_native_inventory(roots):
    evidence = compare_trees(*roots)
    assert len(evidence) == 2
    assert all(row["frozen"]["sha256"] == row["launcher"]["sha256"] for row in evidence)


@pytest.mark.parametrize("name", ["core.py", "kernel.so"])
def test_changed_source_or_library_rejected(roots, name):
    (roots[1] / "orthohmm" / name).write_bytes(b"changed")
    with pytest.raises(ValueError, match="file differs"):
        compare_trees(*roots)


@pytest.mark.parametrize("name", ["extra.py", "extra.so"])
def test_unmanifested_source_or_library_rejected(roots, name):
    (roots[1] / "orthohmm" / name).write_bytes(b"extra")
    with pytest.raises(ValueError, match="file sets differ"):
        compare_trees(*roots)


def test_missing_native_library_rejected(roots):
    (roots[1] / "orthohmm/kernel.so").unlink()
    with pytest.raises(ValueError, match="file sets differ"):
        compare_trees(*roots)

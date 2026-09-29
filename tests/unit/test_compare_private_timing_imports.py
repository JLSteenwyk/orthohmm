import pytest

from benchmark_tools.compare_private_timing_imports import compare, record, run


def test_import_comparison_preserves_differences(tmp_path):
    old_root, new_root = tmp_path / "old/site-packages", tmp_path / "new/site-packages"
    old_root.mkdir(parents=True)
    new_root.mkdir(parents=True)
    for root in (old_root, new_root):
        (root / "same.py").write_text("same")
    (old_root / "changed.py").write_text("before")
    (new_root / "changed.py").write_text("after")
    old = dict(modules={"same": str(old_root / "same.py"), "changed": str(old_root / "changed.py"),
                        "omitted": str(old_root / "hook.py"), "stdlib": "/usr/lib/python/os.py"},
               files=[record(old_root / name) for name in ("same.py", "changed.py")])
    new = dict(modules={"same": str(new_root / "same.py"), "changed": str(new_root / "changed.py"),
                        "added": str(new_root / "additional.py")})
    result = compare(old, new, tmp_path / "core")
    assert result["compared"] == 2 and result["identical"] == 1
    assert [r["module"] for r in result["changed"]] == ["changed"]
    assert [r["module"] for r in result["omitted"]] == ["omitted"]
    assert result["additional_module_names"] == ["added"]


def test_unpinned_import_rejected(tmp_path):
    path = str(tmp_path / "site-packages/unknown.py")
    with pytest.raises(ValueError, match="lacks a hash"):
        compare(dict(modules={"unknown": path}, files=[]), dict(modules={"unknown": path}), tmp_path / "core")


def test_no_overwrite(tmp_path):
    with pytest.raises(FileExistsError):
        run(tmp_path / "missing", "0" * 64, tmp_path / "missing", tmp_path)


def test_wrong_prior_pin(tmp_path):
    path = tmp_path / "prior.json"
    path.write_text("{}")
    with pytest.raises(ValueError, match="checksum"):
        run(path, "0" * 64, tmp_path / "missing", tmp_path / "output")

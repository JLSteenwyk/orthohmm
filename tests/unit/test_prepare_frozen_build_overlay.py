from pathlib import Path
import subprocess

import pytest

from benchmark_tools import prepare_frozen_build_overlay as module

ROOT = Path(__file__).resolve().parents[2]


def test_real_git_staging_changes_only_setup(tmp_path):
    output = tmp_path / "build"
    report = module.prepare(ROOT, output)
    assert len(report["files"]) == 43
    assert sum(r["setup_overlay"] for r in report["files"]) == 1
    for row in report["files"]:
        path = output / "source" / row["relative_path"]
        expected = subprocess.check_output(["git", "-C", str(ROOT), "show",
            (module.BUILD_REVISION if row["setup_overlay"] else module.REVISION) + ":" + row["relative_path"]])
        assert path.read_bytes() == expected
        assert path.stat().st_mode & 0o777 == 0o644
    assert (output / "historical_setup.py").read_bytes() == subprocess.check_output(
        ["git", "-C", str(ROOT), "show", module.REVISION + ":setup.py"])
    assert report["build_executed"] is False
    with pytest.raises(FileExistsError):
        module.prepare(ROOT, output)


@pytest.mark.parametrize("bad", ["missing_setup", "missing_package", "traversal", "absolute"])
def test_invalid_inventory_refused_before_output(tmp_path, monkeypatch, bad):
    item = dict(content=b"", mode=0o644, git_blob="fixture")
    files = {"setup.py": item, "orthohmm/__init__.py": item}
    if bad == "missing_setup":
        files.pop("setup.py")
    elif bad == "missing_package":
        files.pop("orthohmm/__init__.py")
    else:
        files["../escaped.py" if bad == "traversal" else "/escaped.py"] = item
    monkeypatch.setattr(module, "git_files", lambda repo: files)
    with pytest.raises(ValueError):
        module.prepare(ROOT, tmp_path / "out")
    assert not (tmp_path / "out").exists()

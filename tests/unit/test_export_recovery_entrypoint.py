import json

import pytest

from benchmark_tools import export_recovery_entrypoint as module


@pytest.fixture
def exported(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output",
                        lambda command, cwd: command[-1].encode())
    directory = tmp_path / "export"
    module.export(tmp_path, directory)
    return directory


def test_exact_inventory(exported):
    manifest = module.verify(exported)
    assert tuple(manifest["files"]) == module.FILES
    assert manifest["publication_ready"] is False


def test_refuses_existing(exported):
    with pytest.raises(FileExistsError):
        module.export(exported, exported)


def test_changed_file(exported):
    (exported / module.FILES[0]).write_text("changed")
    with pytest.raises(ValueError, match="Changed export"):
        module.verify(exported)


def test_changed_inventory(exported):
    path = exported / "manifest.json"
    manifest = json.loads(path.read_text())
    manifest["files"]["../outside"] = manifest["files"].pop(module.FILES[0])
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="inventory"):
        module.verify(exported)


def test_symlink_rejected(exported, tmp_path):
    path = exported / module.FILES[0]
    original = tmp_path / "original"
    path.rename(original)
    path.symlink_to(original)
    with pytest.raises(ValueError, match="Changed export"):
        module.verify(exported)

import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import export_frozen_figure_helper as helper


@pytest.fixture
def repository(tmp_path, monkeypatch):
    repo = tmp_path / "repository"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    path = repo / helper.SOURCE
    path.parent.mkdir()
    payload = b"# Historical fixture, not current checkout bytes\n"
    path.write_bytes(payload)
    subprocess.run(["git", "add", helper.SOURCE], cwd=repo, check=True)
    subprocess.run(["git", "-c", "user.name=Test", "-c", "user.email=test@example.invalid",
                    "-c", "commit.gpgsign=false", "commit", "-qm", "fixture"], cwd=repo, check=True)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    blob = subprocess.check_output(["git", "rev-parse", f"{commit}:{helper.SOURCE}"], cwd=repo, text=True).strip()
    for key, value in {"COMMIT": commit, "BLOB": blob, "SIZE": len(payload),
                       "SHA256": hashlib.sha256(payload).hexdigest()}.items():
        monkeypatch.setattr(helper, key, value)
    path.write_bytes(b"# Different dirty checkout\n")
    return repo, payload


def test_export_uses_historical_blob_and_relative_mapping(repository, tmp_path):
    repo, payload = repository
    output = tmp_path / "export"
    result = helper.export(repo, output)
    assert (output / helper.SOURCE).read_bytes() == payload
    assert (repo / helper.SOURCE).read_bytes() == b"# Different dirty checkout\n"
    assert json.loads((output / "manifest.json").read_text()) == result
    assert not Path(result["export"]["path"]).is_absolute()
    assert result["publication_ready"] is False
    assert result["historical_repository_relative_path"] == helper.HISTORICAL_PATH


@pytest.mark.parametrize("key,value", [("SHA256", "0" * 64), ("SIZE", 1), ("BLOB", "0" * 40)])
def test_identity_mismatch_leaves_no_output(repository, tmp_path, monkeypatch, key, value):
    repo, _ = repository
    monkeypatch.setattr(helper, key, value)
    output = tmp_path / "export"
    with pytest.raises(ValueError, match="identity mismatch"):
        helper.export(repo, output)
    assert not output.exists()


def test_missing_commit_leaves_no_output(repository, tmp_path, monkeypatch):
    repo, _ = repository
    monkeypatch.setattr(helper, "COMMIT", "0" * 40)
    output = tmp_path / "export"
    with pytest.raises(subprocess.CalledProcessError):
        helper.export(repo, output)
    assert not output.exists()


@pytest.mark.parametrize("kind", ["directory", "file", "dangling_symlink"])
def test_existing_destination_refused(repository, tmp_path, kind):
    repo, _ = repository
    output = tmp_path / "export"
    if kind == "directory":
        output.mkdir()
    elif kind == "file":
        output.write_text("retain me")
    else:
        output.symlink_to(tmp_path / "absent")
    with pytest.raises(FileExistsError):
        helper.export(repo, output)


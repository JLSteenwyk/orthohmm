import json
from pathlib import Path
import shutil
import subprocess

import pytest

from benchmark_tools import bundle_publication_figures as bundle


def commit(repo):
    subprocess.run(["git", "add", "."], cwd=repo, check=True)
    subprocess.run(["git", "-c", "user.name=Test", "-c", "user.email=test@example.invalid",
                    "-c", "commit.gpgsign=false", "commit", "-qm", "fixture"], cwd=repo, check=True)


@pytest.fixture
def repository(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    panel = "synthetic"
    figure_dir = repo / "benchmark_tools/results" / panel
    figure_dir.mkdir(parents=True)
    rows = []
    for extension in ("png", "pdf", "svg"):
        path = figure_dir / ("figure." + extension)
        path.write_bytes(b"test placeholder; not a rendered figure")
        rows.append({"path": str(path), **bundle.identity(path.read_bytes())})
    manifest_path = figure_dir / "manifest.json"
    manifest_path.write_text(json.dumps({"outputs": rows}))
    audit = {"status": "retained_figure_bytes_verified", "panels": [{
        "panel": panel, "status": "all_recorded_bytes_match", "output_count": 3,
        "manifest": {"path": str(manifest_path), **bundle.identity(manifest_path.read_bytes())},
        "files": [{"expected": r, "actual": r, "status": "matches",
                   "repository_relative_path": str(Path(r["path"]).relative_to(repo))} for r in rows],
    }]}
    (repo / bundle.AUDIT).write_text(json.dumps(audit))
    (repo / "LICENSE.md").write_text("Test fixture license placeholder\n")
    shutil.copyfile(bundle.__file__, repo / bundle.RUNNER)
    commit(repo)
    return repo


def test_committed_export_verifies_after_relocation_without_repo(repository, tmp_path):
    out = tmp_path / "export"
    (repository / "benchmark_tools/results/synthetic/figure.png").write_bytes(b"dirty bytes")
    result = bundle.build(repository, "HEAD", out)
    moved = tmp_path / "moved"
    out.rename(moved)
    repository.rename(tmp_path / "unavailable-original-repo")
    assert bundle.verify(moved) == result
    run = subprocess.run(["python", "-I", str(moved / bundle.RUNNER), "verify", str(moved)],
                         cwd=tmp_path, capture_output=True, text=True, check=True)
    assert json.loads(run.stdout) == result
    assert result["panels"] == 1 and result["outputs"] == 3
    assert result["publication_ready"] is False


@pytest.mark.parametrize("mutation", ["bytes", "missing", "extra", "symlink", "mapping", "traversal", "duplicate", "scope"])
def test_verifier_rejects_damage(repository, tmp_path, mutation):
    out = tmp_path / "export"
    bundle.build(repository, "HEAD", out)
    figure = out / "benchmark_tools/results/synthetic/figure.png"
    manifest = json.loads((out / "bundle.json").read_text())
    if mutation == "bytes":
        figure.write_bytes(b"changed")
    elif mutation == "missing":
        figure.unlink()
    elif mutation == "extra":
        (out / "extra").write_bytes(b"unexpected")
    elif mutation == "symlink":
        figure.unlink()
        figure.symlink_to(repository / "benchmark_tools/results/synthetic/figure.png")
    elif mutation == "mapping":
        manifest["panels"][0]["relocations"].popitem()
    elif mutation == "traversal":
        manifest["files"][0]["path"] = "../escape"
    elif mutation == "duplicate":
        manifest["files"].append(manifest["files"][0])
    elif mutation == "scope":
        manifest["publication_ready"] = True
    (out / "bundle.json").write_text(json.dumps(manifest))
    with pytest.raises((ValueError, FileNotFoundError)):
        bundle.verify(out)


def test_committed_changed_evidence_rejected_before_export(repository, tmp_path):
    (repository / "benchmark_tools/results/synthetic/figure.png").write_bytes(b"changed")
    commit(repository)
    out = tmp_path / "export"
    with pytest.raises(ValueError, match="Committed bytes differ"):
        bundle.build(repository, "HEAD", out)
    assert not out.exists()


def test_existing_destination_refused(repository, tmp_path):
    out = tmp_path / "export"
    out.mkdir()
    with pytest.raises(FileExistsError):
        bundle.build(repository, "HEAD", out)


@pytest.mark.parametrize("name", ["../x", "/x", "a/../x", "a//x", "./x", ""])
def test_unsafe_paths(name):
    with pytest.raises(ValueError):
        bundle.safe_relative(name)

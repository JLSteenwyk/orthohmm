from pathlib import Path
import subprocess

import pytest

from benchmark_tools import audit_libleiden_version as audit


@pytest.fixture
def repository(tmp_path):
    def git(*args):
        return subprocess.check_output(["git", "-C", str(tmp_path), *args], stderr=subprocess.STDOUT).decode().strip()
    git("init")
    git("config", "user.name", "Synthetic test")
    git("config", "user.email", "test@example.invalid")
    git("config", "commit.gpgsign", "false")
    git("config", "tag.gpgsign", "false")
    git("remote", "add", "origin", audit.ORIGIN)
    for name in audit.FILES:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("synthetic test source\n")
    (tmp_path / "etc/cmake/version.cmake").write_text(
        'git_describe(PACKAGE_VERSION)\nstring(REGEX MATCH "^[^-]+" PACKAGE_VERSION_BASE "${PACKAGE_VERSION}")\n')
    git("add", ".")
    git("commit", "-m", "old annotated release")
    git("tag", "-a", "0.11.1", "-m", "synthetic old release")
    (tmp_path / "new.txt").write_text("new release content")
    git("add", ".")
    git("commit", "-m", "lightweight release")
    git("tag", "0.12.0")
    return tmp_path, git("rev-parse", "HEAD"), git


def test_lightweight_tag_does_not_override_annotated_describe(repository):
    root, commit, _ = repository
    result = audit.inspect(root, commit)
    assert result["git_describe"].startswith("0.11.1-1-g")
    assert result["git_describe_with_tags"] == "0.12.0"
    assert result["inferred_package_version_base"] == "0.11.1"
    assert not result["wheel_source_identity_established"] and not result["cmake_executed"]
    assert len(result["files"]) == len(audit.FILES) + 1


def test_annotated_release_changes_lookup(repository):
    root, commit, git = repository
    git("tag", "-d", "0.12.0")
    git("tag", "-a", "0.12.0", "-m", "synthetic annotated replacement")
    result = audit.inspect(root, commit)
    assert result["git_describe"] == "0.12.0"
    assert result["inferred_package_version_base"] == "0.12.0"


def test_reads_committed_objects_not_worktree(repository):
    root, commit, _ = repository
    (root / "etc/cmake/version.cmake").write_text("uncommitted corruption")
    assert audit.inspect(root, commit)["inferred_package_version_base"] == "0.11.1"


@pytest.mark.parametrize("kind", ["wrong_origin", "wrong_commit", "override", "symlink", "changed_derivation"])
def test_inconsistent_source_rejected(repository, kind):
    root, commit, git = repository
    if kind == "wrong_origin":
        git("remote", "set-url", "origin", "https://example.invalid/other.git")
    elif kind == "wrong_commit":
        commit = git("rev-parse", "HEAD^")
    else:
        if kind == "override":
            (root / "VERSION").write_text("0.12.0")
        elif kind == "symlink":
            (root / "linked").symlink_to("LICENSE")
        else:
            (root / "etc/cmake/version.cmake").write_text("changed")
        git("add", ".")
        git("commit", "-m", "synthetic bad input")
        git("tag", "-f", "0.12.0")
        commit = git("rev-parse", "HEAD")
    with pytest.raises(ValueError):
        audit.inspect(root, commit)

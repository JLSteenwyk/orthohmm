import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools import bundle_publication_source as module


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


@pytest.fixture
def exported(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    repo.mkdir()
    git(repo, "init", "-q")
    git(repo, "config", "user.name", "Test")
    git(repo, "config", "user.email", "test@example.invalid")
    for name, content in {"LICENSE.md": "Fixture license", "README.md": "Frozen", "requirements.txt": "",
                          "setup.py": "pass\n", "orthohmm/version.py": "VERSION = 'fixture'\n"}.items():
        path = repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    git(repo, "add", ".")
    git(repo, "commit", "-qm", "Scientific")
    scientific = git(repo, "rev-parse", "HEAD")
    original = module.SCIENTIFIC_REVISION
    monkeypatch.setattr(module, "SCIENTIFIC_REVISION", scientific)
    runner = Path(module.__file__).read_text().replace(original, scientific)
    for name, content in {module.RUNNER: runner, module.GUIDE: "Source scope", "benchmark_tools/probe.py": "pass\n",
                          "tests/unit/test_example.py": "def test_example():\n    pass\n",
                          "benchmark_tools/results/private.json": '{"omitted": true}',
                          "tests/samples/raw.fa": ">gene\nA\n", "orthohmm/version.py": "VERSION = 'changed'\n"}.items():
        path = repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    git(repo, "add", ".")
    git(repo, "commit", "-qm", "Workflow")
    revision = git(repo, "rev-parse", "HEAD")
    # Mutable source must never enter the committed export.
    (repo / "benchmark_tools/probe.py").write_text("not valid python !")
    bundle = tmp_path / "bundle"
    result = module.build(repo, revision, bundle)
    return repo, bundle, result


def rewrite_index(bundle, change):
    path = bundle / "SOURCE_INDEX.json"
    value = json.loads(path.read_text())
    change(value)
    path.write_text(json.dumps(value))
    return module.identity(path.read_bytes())["sha256"]


def test_component_revisions_and_exclusions(exported):
    repo, bundle, result = exported
    assert result["status"] == "publication_source_components_verified"
    assert (bundle / "scientific/orthohmm/version.py").read_text() == "VERSION = 'fixture'\n"
    assert (bundle / "workflow/benchmark_tools/probe.py").read_text() == "pass\n"
    assert not (bundle / "workflow/benchmark_tools/results").exists()
    assert not (bundle / "workflow/tests/samples").exists()
    assert result["components"] == dict(scientific=5, workflow=5)
    assert result["executable_benchmark_reproduced"] is False
    assert result["publication_ready"] is False
    assert result["redistribution_clearance"] is False


def test_relocated_verifier_without_checkout_or_git(exported, tmp_path):
    repo, bundle, result = exported
    relocated = tmp_path / "relocated"
    shutil.copytree(bundle, relocated)
    shutil.rmtree(bundle)
    shutil.rmtree(repo)
    command = [sys.executable, "-I", "-B", str(relocated / "workflow" / module.RUNNER),
               "verify", str(relocated), "--manifest-sha256", result["manifest"]["sha256"]]
    checked = subprocess.run(command, cwd=tmp_path, env={"PATH": "/no-git-here"},
                             capture_output=True, text=True, timeout=10)
    assert checked.returncode == 0, checked.stderr
    assert json.loads(checked.stdout) == result


@pytest.mark.parametrize("name", ["../escape", "/absolute", "a/../b", "a//b", "a/./b", "", ".", 1])
def test_bad_paths_rejected(name):
    with pytest.raises(ValueError):
        module.relative(name)


@pytest.mark.parametrize("defect", ["bytes", "missing", "extra", "symlink", "mode", "index_symlink"])
def test_changed_payload_rejected(exported, defect):
    repo, bundle, result = exported
    path = bundle / "workflow/benchmark_tools/probe.py"
    if defect == "bytes":
        path.write_text("pass # changed\n")
    elif defect == "missing":
        path.unlink()
    elif defect == "extra":
        (bundle / "extra.txt").write_text("extra")
    elif defect == "symlink":
        path.unlink()
        path.symlink_to(repo / "setup.py")
    elif defect == "mode":
        path.chmod(0o777)
    else:
        index = bundle / "SOURCE_INDEX.json"
        outside = bundle.parent / "outside.json"
        index.rename(outside)
        index.symlink_to(outside)
    with pytest.raises((ValueError, FileNotFoundError)):
        module.verify(bundle, result["manifest"]["sha256"])


@pytest.mark.parametrize("defect", ["duplicate", "escape", "mapping", "revision", "ready", "clearance", "omitted"])
def test_invalid_index_rejected_even_with_new_digest(exported, defect):
    repo, bundle, result = exported
    def change(value):
        if defect == "duplicate":
            value["files"].append(value["files"][0])
        elif defect == "escape":
            value["files"][0]["path"] = "../escape"
        elif defect == "mapping":
            value["files"][0]["git_path"] = "not-the-file"
        elif defect == "revision":
            value["scientific_revision"] = "f" * 40
        elif defect == "ready":
            value["publication_ready"] = True
        elif defect == "clearance":
            value["redistribution_clearance"] = True
        else:
            value["files"].pop()
    digest = rewrite_index(bundle, change)
    with pytest.raises(ValueError):
        module.verify(bundle, digest)


def test_external_digest_rejects_coordinated_change(exported):
    repo, bundle, result = exported
    path = bundle / "workflow/benchmark_tools/probe.py"
    path.write_text("pass # changed\n")
    def change(value):
        for row in value["files"]:
            if row["path"] == "workflow/benchmark_tools/probe.py":
                row.update(module.identity(path.read_bytes()))
    rewrite_index(bundle, change)
    with pytest.raises(ValueError, match="external anchor"):
        module.verify(bundle, result["manifest"]["sha256"])


def test_existing_destination_refused(exported):
    repo, bundle, result = exported
    with pytest.raises(FileExistsError):
        module.build(repo, result["workflow_revision"], bundle)

import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools import bundle_publication_source as module


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


@pytest.fixture(params=["source-only", "orthobench-inputs", "native-preparation", "native-wheels"])
def exported(tmp_path, monkeypatch, request):
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
    if request.param != "source-only":
        root = Path(module.__file__).resolve().parent.parent
        support_names = {*module.support_pins(request.param), "benchmark_tools/verify_orthobench_acquisition.py",
                         "benchmark_tools/rebind_orthobench_data.py"}
        if request.param in {"native-preparation", "native-wheels"}:
            support_names |= module.BASE_HELPERS
        if request.param == "native-wheels":
            support_names |= module.WHEEL_HELPERS
        for name in support_names:
            target = repo / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes((root / name).read_bytes())
    git(repo, "add", ".")
    git(repo, "commit", "-qm", "Workflow")
    revision = git(repo, "rev-parse", "HEAD")
    # Mutable source must never enter the committed export.
    (repo / "benchmark_tools/probe.py").write_text("not valid python !")
    bundle = tmp_path / "bundle"
    result = module.build(repo, revision, bundle, request.param)
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
    index = json.loads((bundle / "SOURCE_INDEX.json").read_text())
    support = index.get("profile") in module.PROFILES - {"source-only"}
    assert (bundle / "workflow/benchmark_tools/results").exists() is support
    if support:
        assert {p.relative_to(bundle / "workflow").as_posix()
                for p in (bundle / "workflow/benchmark_tools/results").iterdir()} == set(module.support_pins(index["profile"]))
    else:
        assert index["schema"] == "publication_source_components_v1"
    assert not (bundle / "workflow/tests/samples").exists()
    workflow_count = 5 if not support else 10
    if index.get("profile") in {"native-preparation", "native-wheels"}:
        workflow_count += 1 + len(module.BASE_HELPERS)
    if index.get("profile") == "native-wheels":
        workflow_count += len(module.WHEEL_SUPPORT_PINS) + len(module.WHEEL_HELPERS - module.BASE_HELPERS)
    assert result["components"] == dict(scientific=5, workflow=workflow_count)
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


def test_unknown_profile_rejected(tmp_path):
    with pytest.raises(ValueError, match="Unknown source profile"):
        module.build(tmp_path, "HEAD", tmp_path / "output", "everything")


def test_default_selection_still_excludes_present_support(exported, tmp_path):
    repo, _, result = exported
    bundle = tmp_path / "default-export"
    module.build(repo, result["workflow_revision"], bundle)
    index = json.loads((bundle / "SOURCE_INDEX.json").read_text())
    assert index["schema"] == "publication_source_components_v1"
    assert "profile" not in index
    assert not (bundle / "workflow/benchmark_tools/results").exists()


@pytest.mark.parametrize("defect", ["profile", "schema", "support_pin", "support_missing", "helper_missing"])
def test_support_contract_rejected(exported, defect):
    _, bundle, result = exported
    path = bundle / "SOURCE_INDEX.json"
    index = json.loads(path.read_text())
    if index.get("profile") not in module.PROFILES - {"source-only"}:
        # v1 cannot reinterpret the legacy selection even with a new external digest.
        digest = rewrite_index(bundle, lambda value: value.update(profile="orthobench-inputs"))
        with pytest.raises(ValueError, match="Historical source schema"):
            module.verify(bundle, digest)
        return
    if defect in {"profile", "schema"}:
        index["profile" if defect == "profile" else "schema"] = "unreviewed"
    else:
        name = ("benchmark_tools/rebind_orthobench_data.py" if defect == "helper_missing"
                else next(iter(module.SUPPORT_PINS)))
        target = bundle / "workflow" / name
        if defect == "support_pin":
            target.write_text("{}\n")
            for row in index["files"]:
                if row["path"] == "workflow/" + name:
                    row.update(module.identity(target.read_bytes()))
        else:
            target.unlink()
            index["files"] = [row for row in index["files"] if row["path"] != "workflow/" + name]
    path.write_text(json.dumps(index))
    with pytest.raises(ValueError):
        module.verify(bundle, module.identity(path.read_bytes())["sha256"])


@pytest.mark.parametrize("defect", ["missing", "changed"])
def test_builder_rejects_changed_support_before_output(exported, tmp_path, defect):
    repo, bundle, result = exported
    name = next(iter(module.SUPPORT_PINS))
    target = repo / name
    if target.exists() and defect == "missing":
        git(repo, "rm", name)
    elif defect == "changed":
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("{}\n")
        git(repo, "add", name)
    else:
        # The legacy fixture already lacks the acquisition manifests.
        pass
    if git(repo, "diff", "--cached", "--name-only"):
        git(repo, "commit", "-qm", "Defective support")
    output = tmp_path / "rejected"
    with pytest.raises(ValueError):
        module.build(repo, "HEAD", output, "orthobench-inputs")
    assert not output.exists()


def test_native_base_manifest_pin_rejected(exported):
    _, bundle, _ = exported
    path = bundle / "SOURCE_INDEX.json"
    index = json.loads(path.read_text())
    if index.get("profile") not in {"native-preparation", "native-wheels"}:
        assert not (bundle / "workflow" / next(iter(module.BASE_SUPPORT_PINS))).exists()
        return
    name = next(iter(module.BASE_SUPPORT_PINS))
    target = bundle / "workflow" / name
    target.write_text("{}\n")
    for row in index["files"]:
        if row["path"] == "workflow/" + name:
            row.update(module.identity(target.read_bytes()))
    path.write_text(json.dumps(index))
    with pytest.raises(ValueError, match="Frozen acquisition/runtime"):
        module.verify(bundle, module.identity(path.read_bytes())["sha256"])


@pytest.mark.parametrize("name", sorted(module.WHEEL_SUPPORT_PINS))
def test_native_wheel_support_pin_rejected(exported, name):
    _, bundle, _ = exported
    index_path = bundle / "SOURCE_INDEX.json"
    index = json.loads(index_path.read_bytes())
    target = bundle / "workflow" / name
    if index.get("profile") != "native-wheels":
        assert not target.exists()
        return
    target.write_text("changed\n")
    for row in index["files"]:
        if row["path"] == "workflow/" + name:
            row.update(module.identity(target.read_bytes()))
    index_path.write_text(json.dumps(index))
    with pytest.raises(ValueError, match="Frozen acquisition/runtime"):
        module.verify(bundle, module.identity(index_path.read_bytes())["sha256"])

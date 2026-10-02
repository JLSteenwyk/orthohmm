import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

import pytest

from benchmark_tools import bundle_swiss_descriptive_tables as bundle

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def source_repo(tmp_path):
    repo = tmp_path / "source repo with spaces"
    sources = [bundle.RUNNER, bundle.CHECKER, "LICENSE.md",
               *["benchmark_tools/results/" + name for name in bundle.data_paths()]]
    for name in sources:
        target = repo / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    subprocess.run(["git", "init", "--quiet", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "add", "--", *sources], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid",
                    "commit", "--quiet", "-m", "Fixture component source"], check=True)
    return repo


@pytest.fixture
def component(source_repo, tmp_path):
    output = tmp_path / "review component with spaces"
    result = bundle.build(source_repo.resolve(), "HEAD", output.resolve())
    return output, result["manifest"]["sha256"], result


@pytest.mark.parametrize("name", ["../escape", "/absolute", "results//x", "results/../x", "C:/x",
                                  "results\\evil", "results/./x"])
def test_unsafe_member_names_rejected(name):
    with pytest.raises(ValueError, match="Unsafe component member"):
        bundle.safe_name(name)


def test_committed_export_and_reproduction(component, tmp_path):
    directory, digest, built = component
    checked, payloads = bundle.verify(directory, digest)
    assert checked == built and checked["files"] == len(payloads) == 21
    assert checked["numerical_reproduction_executed"] is checked["publication_ready"] is False
    output = (tmp_path / "result.json").resolve()
    result = bundle.reproduce(directory, digest, output)
    assert result["rows"] == 208 and result["score_and_difference_cells"] == 984
    assert result["component"] == checked and result["raw_inputs_revalidated"] is False
    assert result["component_numerical_reproduction_executed"] is True
    assert json.loads(output.read_text())["checker"]["sha256"] == bundle.CHECKER_SHA
    original = output.read_bytes()
    with pytest.raises(FileExistsError):
        bundle.reproduce(directory, digest, output)
    assert output.read_bytes() == original
    with pytest.raises(ValueError, match="must not mutate"):
        bundle.reproduce(directory, digest, directory / "new.json")


def test_export_ignores_dirty_source_payload(source_repo, tmp_path):
    path = source_repo / "benchmark_tools/results/qfo_recovered_swiss_uncertainty_22178.json"
    path.write_bytes(b"not the committed payload")
    result = bundle.build(source_repo.resolve(), "HEAD", (tmp_path / "component").resolve())
    assert result["files"] == 21
    with pytest.raises(FileExistsError):
        bundle.build(source_repo.resolve(), "HEAD", (tmp_path / "component").resolve())


def test_changed_committed_builder_rejected_before_export(source_repo, tmp_path):
    path = source_repo / bundle.RUNNER
    path.write_bytes(path.read_bytes() + b"\n# Different committed builder\n")
    subprocess.run(["git", "-C", str(source_repo), "add", bundle.RUNNER], check=True)
    subprocess.run(["git", "-C", str(source_repo), "-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid",
                    "commit", "--quiet", "-m", "Changed fixture builder"], check=True)
    output = (tmp_path / "component").resolve()
    with pytest.raises(ValueError, match="Executing builder"):
        bundle.build(source_repo.resolve(), "HEAD", output)
    assert not output.exists()


def test_archive_restores_full_arithmetic(component, tmp_path):
    directory, digest, built = component
    archive = (tmp_path / "tables.tar.gz").resolve()
    result = bundle.archive(directory, digest, archive)
    assert result["members"] == 22 and result["component"] == built
    restored = bundle.restore_and_reproduce(archive, digest, (tmp_path / "restored.json").resolve())
    assert restored["rows"] == 208 and restored["score_and_difference_cells"] == 984
    assert restored["component"] == built
    assert not Path(restored["checker"]["path"]).exists()
    with pytest.raises(FileExistsError):
        bundle.archive(directory, digest, archive)
    with pytest.raises(ValueError, match="must not mutate"):
        bundle.archive(directory, digest, directory / "tables.tar.gz")


@pytest.mark.parametrize("fault", ["manifest", "digest", "payload", "extra", "missing", "symlink"])
def test_changed_component_is_rejected(component, fault, tmp_path):
    directory, digest, _ = component
    if fault == "manifest":
        p = directory / "bundle.json"
        p.write_bytes(p.read_bytes() + b"\n")
    elif fault == "digest":
        digest = "0" * 64
    elif fault == "payload":
        p = directory / "results/qfo_recovered_swiss_uncertainty_22178.json"
        p.write_bytes(p.read_bytes() + b"\n")
    elif fault == "extra":
        (directory / "extra.txt").write_text("extra")
    elif fault == "missing":
        (directory / "README.md").unlink()
    else:
        p = directory / "LICENSE.md"
        target = tmp_path / "license.txt"
        p.rename(target)
        p.symlink_to(target)
    with pytest.raises((ValueError, FileNotFoundError)):
        bundle.verify(directory, digest)


@pytest.mark.parametrize("key,value", [("scope", "native_inference"), ("publication_ready", True),
                                     ("raw_source_admission", True), ("source_commit", "not-a-commit")])
def test_wrong_scope_rejected_even_with_rehashed_index(component, key, value):
    directory, _, _ = component
    path = directory / "bundle.json"
    manifest = json.loads(path.read_text())
    manifest[key] = value
    path.write_text(json.dumps(manifest))
    digest = bundle.identity(path.read_bytes())["sha256"]
    with pytest.raises(ValueError, match="Wrong component scope"):
        bundle.verify(directory, digest)


@pytest.mark.parametrize("fault", ["traversal", "symlink", "duplicate", "executable", "oversized", "missing"])
def test_unsafe_or_incomplete_archive_rejected_before_arithmetic(tmp_path, fault):
    archive = (tmp_path / "bad.tar.gz").resolve()
    info = tarfile.TarInfo("../escape" if fault == "traversal" else "README.md")
    info.mode, info.size = 0o644, 1
    if fault == "symlink":
        info.type, info.linkname, info.size = tarfile.SYMTYPE, "/outside", 0
    elif fault == "executable":
        info.mode = 0o755
    elif fault == "oversized":
        info.size = bundle.MAX_BYTES + 1
    with tarfile.open(archive, "w:gz") as handle:
        handle.addfile(info, None if fault == "oversized" else io.BytesIO(b"x"))
        if fault == "duplicate":
            handle.addfile(info, io.BytesIO(b"x"))
    output = (tmp_path / "result.json").resolve()
    with pytest.raises(ValueError):
        bundle.restore_and_reproduce(archive, "0"*64, output)
    assert not output.exists() and not (tmp_path / "escape").exists()

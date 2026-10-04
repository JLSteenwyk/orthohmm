"""Local versioned package integrity/relocation; no scientific jobs execute."""

import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

import pytest

from benchmark_tools import bundle_publication_package as package


def write(path, data):
    path.write_text(json.dumps(data))
    return package.identity(path)["sha256"]


@pytest.fixture
def selected(tmp_path):
    root = tmp_path / "inputs"
    root.mkdir()
    (root / "README.md").write_text("Local working candidate, not publication readiness.")
    (root / "payload.bin").write_bytes(b"retained artifact" * 180000)
    plan = {"schema": "publication_package_selection_v1", "version": "orthohmm-study-2026.10.04-rc1",
            "scientific_revision": "a" * 40, "workflow_revision": "b" * 40,
            "publication_ready": False, "public_archive_uploaded": False,
            "limitations": ["Synthetic fixture, not scientific evidence."],
            "files": [{"target": name, "source": name, **package.identity(root / name)}
                      for name in ("README.md", "payload.bin")]}
    selection = root / "selection.json"
    anchor = write(selection, plan)
    return root, selection, anchor, plan


@pytest.mark.parametrize("optimized", [False, True])
def test_archive_restore_and_copied_cli(selected, tmp_path, optimized):
    root, selection, anchor, _ = selected
    built = package.build(root, selection, anchor, tmp_path / "built")
    assert built["files"] == 4
    archive = package.archive(tmp_path / "built", built["manifest"]["sha256"], tmp_path / "candidate.tar.gz")
    restored = package.restore(tmp_path / "candidate.tar.gz", archive["archive"]["sha256"],
                               built["manifest"]["sha256"], tmp_path / "restored")
    assert restored == built
    completed = subprocess.run([sys.executable, "-I", "-S", "-B", *(["-O"] if optimized else []),
        str(tmp_path / "restored" / package.RUNNER), "verify", str(tmp_path / "restored"),
        "--manifest-sha256", built["manifest"]["sha256"]], cwd=tmp_path, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout) == built
    assert completed.stderr == ""
    again = package.archive(tmp_path / "built", built["manifest"]["sha256"], tmp_path / "again.tar.gz")
    assert again["archive"]["sha256"] == archive["archive"]["sha256"]


@pytest.mark.parametrize("damage", ["selection", "payload", "unsafe_source", "unsafe_target", "duplicate", "reserved", "collision", "publication_claim", "revision", "version"])
def test_build_rejects_bad_inputs_before_output(selected, tmp_path, damage):
    root, selection, anchor, plan = selected
    if damage == "selection":
        selection.write_text(selection.read_text() + " ")
    elif damage == "payload":
        (root / "payload.bin").write_bytes(b"changed")
    else:
        if damage == "unsafe_source":
            plan["files"][1]["source"] = "../payload.bin"
        elif damage == "unsafe_target":
            plan["files"][1]["target"] = "../payload.bin"
        elif damage == "duplicate":
            plan["files"][1]["target"] = "README.md"
        elif damage == "reserved":
            plan["files"][1]["target"] = package.INDEX
        elif damage == "collision":
            plan["files"][1]["target"] = "README.md/child"
        elif damage == "publication_claim":
            plan["publication_ready"] = True
        elif damage == "revision":
            plan["workflow_revision"] = "HEAD"
        else:
            plan["version"] = "latest"
        anchor = write(selection, plan)
    with pytest.raises(ValueError):
        package.build(root, selection, anchor, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()


@pytest.mark.parametrize("damage", ["anchor", "missing", "changed", "extra", "mode", "symlink"])
def test_verify_refuses_changed_inventory(selected, tmp_path, damage):
    root, selection, anchor, _ = selected
    built = package.build(root, selection, anchor, tmp_path / "built")
    digest = built["manifest"]["sha256"]
    target = tmp_path / "built/payload.bin"
    if damage == "anchor":
        digest = "0" * 64
    elif damage == "missing":
        target.unlink()
    elif damage == "changed":
        target.write_bytes(b"changed")
    elif damage == "extra":
        (tmp_path / "built/extra").write_text("extra")
    elif damage == "mode":
        target.chmod(0o600)
    else:
        target.unlink()
        target.symlink_to(root / "payload.bin")
    with pytest.raises(ValueError):
        package.verify(tmp_path / "built", digest)


@pytest.mark.parametrize("stage", ["build", "archive", "restore"])
def test_existing_output_preserved(selected, tmp_path, stage):
    root, selection, anchor, _ = selected
    built = package.build(root, selection, anchor, tmp_path / "built")
    arc = package.archive(tmp_path / "built", built["manifest"]["sha256"], tmp_path / "candidate.tar.gz")
    output = tmp_path / "existing"
    output.mkdir()
    marker = output / "marker"
    marker.write_text("preserve")
    with pytest.raises(FileExistsError):
        if stage == "build":
            package.build(root, selection, anchor, output)
        elif stage == "archive":
            package.archive(tmp_path / "built", built["manifest"]["sha256"], output)
        else:
            package.restore(tmp_path / "candidate.tar.gz", arc["archive"]["sha256"], built["manifest"]["sha256"], output)
    assert list(output.iterdir()) == [marker]
    assert marker.read_text() == "preserve"


@pytest.mark.parametrize("damage", ["anchor", "changed", "traversal", "symlink", "duplicate"])
def test_archive_damage_refused_before_extraction(selected, tmp_path, damage):
    root, selection, anchor, _ = selected
    built = package.build(root, selection, anchor, tmp_path / "built")
    arc = package.archive(tmp_path / "built", built["manifest"]["sha256"], tmp_path / "candidate.tar.gz")
    archive_path = tmp_path / "candidate.tar.gz"
    archive_anchor = arc["archive"]["sha256"]
    if damage == "anchor":
        archive_anchor = "0" * 64
    else:
        rewritten = tmp_path / "damaged.tar.gz"
        with tarfile.open(archive_path, "r:gz") as source, tarfile.open(rewritten, "w:gz") as destination:
            for member in source.getmembers():
                content = source.extractfile(member).read()
                if member.name == "payload.bin":
                    if damage == "changed":
                        content = b"changed"
                        member.size = len(content)
                    elif damage == "traversal":
                        member.name = "../payload.bin"
                    elif damage == "symlink":
                        member.type, member.linkname, member.size = tarfile.SYMTYPE, "../outside", 0
                    elif damage == "duplicate":
                        destination.addfile(member, io.BytesIO(content))
                destination.addfile(member, io.BytesIO(content) if member.isfile() else None)
        archive_path = rewritten
        archive_anchor = package.identity(rewritten)["sha256"]
    with pytest.raises(ValueError):
        package.restore(archive_path, archive_anchor, built["manifest"]["sha256"], tmp_path / "bad")
    assert not (tmp_path / "bad").exists()

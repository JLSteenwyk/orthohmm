from io import BytesIO
import hashlib
import json
import shutil
import tarfile

import pytest

from benchmark_tools import export_bundled_source_notices as notices


def fixture(tmp_path, defect=None):
    path = tmp_path / "source.tar.xz"
    rows = []
    with tarfile.open(path, "w:xz") as archive:
        for i, name in enumerate(notices.LLVM_SELECTION):
            data = f"synthetic notice {i}\n".encode()
            rows.append(dict(member=name, bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
            if i == 0 and defect == "missing":
                continue
            entry = tarfile.TarInfo(name)
            entry.size = len(data)
            if i == 0 and defect == "symlink":
                entry.type, entry.linkname, entry.size = tarfile.SYMTYPE, "/outside", 0
            archive.addfile(entry, BytesIO(data))
            if i == 0 and defect == "duplicate":
                archive.addfile(entry, BytesIO(data))
    identity = notices.record(path)
    recipe = tmp_path / "recipe.json"
    recipe.write_text(json.dumps(dict(rendered_source=[
        dict(sha256=identity["sha256"], url="https://example.test/source")])))
    inventory = dict(status="recipe_bound_llvm_source_acquired", source_archive=identity,
                     recipe_receipt=notices.record(recipe), source_url="https://example.test/source",
                     notice_candidates=rows)
    source = tmp_path / "inventory.json"
    source.write_text(json.dumps(inventory))
    return source


def test_llvm_export_verifies_after_originals_removed(tmp_path):
    source = fixture(tmp_path)
    output = tmp_path / "export"
    result = notices.export_llvm(source, output)
    assert result["files"] == 3
    assert result["redistribution_clearance"] is result["publication_ready"] is False
    relocated = tmp_path / "relocated"
    shutil.copytree(output, relocated)
    for name in ("source.tar.xz", "inventory.json", "recipe.json"):
        (tmp_path / name).unlink()
    assert notices.verify(relocated, result["index"]["sha256"])["files"] == 3
    with pytest.raises(FileExistsError):
        notices.export_llvm(source, output)


@pytest.mark.parametrize("defect", ["missing", "duplicate", "symlink"])
def test_invalid_archive_cannot_finalize(tmp_path, defect):
    with pytest.raises(ValueError):
        notices.export_llvm(fixture(tmp_path, defect), tmp_path / "export")
    assert not (tmp_path / "export/SOURCE_NOTICE_INDEX.json").exists()


@pytest.mark.parametrize("defect", ["recipe", "archive", "notice", "duplicate", "omitted", "url"])
def test_changed_binding_rejected(tmp_path, defect):
    source = fixture(tmp_path)
    data = json.loads(source.read_text())
    if defect in ("recipe", "archive"):
        path = tmp_path / ("recipe.json" if defect == "recipe" else "source.tar.xz")
        with path.open("ab") as handle:
            handle.write(b"changed")
    elif defect == "notice":
        data["notice_candidates"][0]["sha256"] = "0" * 64
    elif defect == "duplicate":
        data["notice_candidates"].append(data["notice_candidates"][0])
    elif defect == "omitted":
        data["notice_candidates"].pop()
    else:
        data["source_url"] = "https://example.test/different"
    source.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        notices.export_llvm(source, tmp_path / "export")
    assert not (tmp_path / "export/SOURCE_NOTICE_INDEX.json").exists()


@pytest.mark.parametrize("defect", ["payload", "extra", "symlink", "index"])
def test_relocated_tampering_rejected(tmp_path, defect):
    output = tmp_path / "export"
    result = notices.export_llvm(fixture(tmp_path), output)
    target = output / "llvm/LICENSE.TXT"
    if defect == "payload":
        target.write_text("changed")
    elif defect == "extra":
        (output / "extra").write_text("unlisted")
    elif defect == "index":
        (output / "SOURCE_NOTICE_INDEX.json").write_text("{}")
    else:
        copied = tmp_path / "copy"
        copied.write_bytes(target.read_bytes())
        target.unlink()
        target.symlink_to(copied)
    with pytest.raises(ValueError):
        notices.verify(output, result["index"]["sha256"])

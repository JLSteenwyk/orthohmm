from io import BytesIO
import json
import shutil
import tarfile

import pytest

from benchmark_tools import export_bundled_source_notices as notices


def fixture(tmp_path, *, missing=False, duplicate=False, symlink=False):
    archive = tmp_path / "source.tar.gz"
    selected = []
    with tarfile.open(archive, "w:gz") as tar:
        for i, name in enumerate(notices.SELECTION):
            payload = ("synthetic notice " + str(i)).encode()
            path = tmp_path / ("notice" + str(i))
            path.write_bytes(payload)
            selected.append(dict(member=name, file=notices.record(path)))
            if missing and i == 0:
                continue
            member = tarfile.TarInfo(name)
            member.size = len(payload)
            if symlink and i == 0:
                member.type = tarfile.SYMTYPE
                member.linkname = "/outside"
                member.size = 0
            tar.addfile(member, BytesIO(payload))
            if duplicate and i == 0:
                tar.addfile(member, BytesIO(payload))
    ref = notices.record(archive)
    support = tmp_path / "support.json"
    support.write_text('{"synthetic_only": true}')
    support_ref = notices.record(support)
    inventory = dict(status="distribution_source_packages_acquired_and_inventoried",
        signature_receipt=support_ref, source_archive_inventory=support_ref,
        source_packages=[dict(download=dict(file=ref), members=[dict(member=archive.name, file=ref)])],
        selected_source_material=[dict(archive=ref, selected_members=selected)])
    path = tmp_path / "inventory.json"
    path.write_text(json.dumps(inventory))
    return path


def test_export_and_relocated_offline_verification(tmp_path):
    source = fixture(tmp_path)
    output = tmp_path / "export"
    receipt = notices.export(source, output)
    assert receipt["files"] == 4 and not receipt["redistribution_clearance"]
    moved = tmp_path / "moved"
    shutil.copytree(output, moved)
    source.unlink()
    (tmp_path / "source.tar.gz").unlink()
    other = notices.verify(moved, receipt["index"]["sha256"])
    assert other["index"]["sha256"] == receipt["index"]["sha256"]
    with pytest.raises(FileExistsError):
        notices.export(source, output)


@pytest.mark.parametrize("option", ["missing", "duplicate", "symlink"])
def test_bad_archive_member_rejected(tmp_path, option):
    source = fixture(tmp_path, **{option: True})
    with pytest.raises(ValueError):
        notices.export(source, tmp_path / "export")
    assert not (tmp_path / "export/SOURCE_NOTICE_INDEX.json").exists()


@pytest.mark.parametrize("change", ["payload", "extra", "index", "symlink", "directory_link"])
def test_tampered_relocated_export_rejected(tmp_path, change):
    source = fixture(tmp_path)
    output = tmp_path / "export"
    receipt = notices.export(source, output)
    path = output / "libxml2/Copyright"
    if change == "payload": path.write_text("changed")
    elif change == "extra": (output / "extra").write_text("unlisted")
    elif change == "index": (output / "SOURCE_NOTICE_INDEX.json").write_text("{}")
    elif change == "symlink":
        payload = tmp_path / "copy"
        payload.write_bytes(path.read_bytes())
        path.unlink()
        path.symlink_to(payload)
    else:
        moved = tmp_path / "libxml2"
        (output / "libxml2").rename(moved)
        (output / "libxml2").symlink_to(moved)
    with pytest.raises(ValueError):
        notices.verify(output, receipt["index"]["sha256"])


@pytest.mark.parametrize("change", ["archive", "support", "expected", "unbound", "omitted", "duplicate_record"])
def test_changed_source_evidence_rejected(tmp_path, change):
    source = fixture(tmp_path)
    data = json.loads(source.read_text())
    if change in {"archive", "support"}:
        path = tmp_path / ("source.tar.gz" if change == "archive" else "support.json")
        with path.open("ab") as f: f.write(b"changed")
    else:
        group = data["selected_source_material"][0]
        if change == "expected": group["selected_members"][0]["file"]["sha256"] = "0" * 64
        elif change == "unbound": data["source_packages"][0]["members"] = []
        elif change == "omitted": group["selected_members"].pop()
        else: group["selected_members"].append(group["selected_members"][0])
        source.write_text(json.dumps(data))
    with pytest.raises((ValueError, RuntimeError)):
        notices.export(source, tmp_path / "export")
    assert not (tmp_path / "export/SOURCE_NOTICE_INDEX.json").exists()

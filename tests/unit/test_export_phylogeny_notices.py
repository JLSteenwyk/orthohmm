from io import BytesIO
import json
import os
from pathlib import Path
import shutil
import tarfile

import pytest

from benchmark_tools import export_phylogeny_notices as notices


def fixture(tmp_path, monkeypatch, archive_change=None, header=b"/* synthetic notice */\ncode"):
    sources = tmp_path / "sources"
    sources.mkdir()
    archive = sources / "mafft.tgz"
    refs = []
    with tarfile.open(archive, "w:gz") as stream:
        for name in ("license", "license.extensions"):
            member_name = notices.ROOT + "/" + name
            payload = ("synthetic " + name).encode()
            source = sources / notices.ROOT / name
            source.parent.mkdir(exist_ok=True)
            source.write_bytes(payload)
            refs.append(notices.identity(source))
            if archive_change == "missing" and name == "license":
                continue
            member = tarfile.TarInfo(member_name)
            member.size = len(payload)
            if archive_change == "symlink" and name == "license":
                member.type, member.linkname, member.size = tarfile.SYMTYPE, "/outside", 0
            stream.addfile(member, BytesIO(payload))
            if archive_change == "duplicate" and name == "license":
                stream.addfile(member, BytesIO(payload))
        if archive_change in {"traversal", "absolute", "wrong_root", "oversized", "special"}:
            name = {"traversal": notices.ROOT + "/../bad", "absolute": "/bad",
                    "wrong_root": "another/bad"}.get(archive_change, notices.ROOT + "/bad")
            member = tarfile.TarInfo(name)
            if archive_change == "special":
                member.type = tarfile.FIFOTYPE
            if archive_change == "oversized":
                member.size = 20_000_001
                stream.addfile(member, BytesIO(b"x" * member.size))
            else:
                stream.addfile(member)
    archive_ref = notices.identity(archive)
    monkeypatch.setattr(notices, "MAFFT_SHA", archive_ref["sha256"])
    mafft = sources / "mafft.json"
    mafft.write_text(json.dumps(dict(
        status="core_built_and_installed_phylogeny_fixture_passed", acquisition_url=notices.URL,
        extensions_built=False, archive=archive_ref, source_files=refs)))
    acquired = []
    hashes = dict(notices.FASTTREE_FILES)
    for name, payload in (("LICENSE", b"synthetic separate license"), ("FastTree.c", header)):
        path = sources / name
        path.write_bytes(payload)
        ref = notices.identity(path)
        hashes[name] = ref["sha256"]
        acquired.append(dict(ref, url=notices.FASTTREE_BASE + name))
    monkeypatch.setattr(notices, "FASTTREE_FILES", hashes)
    fasttree = sources / "fasttree.json"
    fasttree.write_text(json.dumps(dict(
        status="source_notices_acquired_and_installed_binary_matched",
        revision=notices.FASTTREE_REVISION, acquired_files=acquired)))
    return [mafft, notices.identity(mafft)["sha256"], fasttree,
            notices.identity(fasttree)["sha256"]]


def rewrite(arguments, index, mutate):
    path = arguments[index]
    data = json.loads(path.read_text())
    mutate(data)
    path.write_text(json.dumps(data))
    arguments[index + 1] = notices.identity(path)["sha256"]


def test_exact_notice_export_and_offline_relocation(tmp_path, monkeypatch):
    arguments = fixture(tmp_path, monkeypatch)
    output = tmp_path / "export"
    receipt = notices.export(*arguments, output)
    assert receipt["files"] == 4 and receipt["redistribution_clearance"] is False
    index = json.loads((output / "SOURCE_NOTICE_INDEX.json").read_bytes())
    header = next(row for row in index["files"] if "leading-comment" in row["member"])
    assert header["byte_range"] == [0, len(b"/* synthetic notice */\n")]
    assert (output / "fasttree/source-header.txt").read_bytes() == b"/* synthetic notice */\n"
    assert (output / "fasttree/LICENSE").read_bytes() == b"synthetic separate license"
    moved = tmp_path / "moved"
    shutil.copytree(output, moved)
    shutil.rmtree(tmp_path / "sources")
    result = notices.verify(moved, receipt["index"]["sha256"])
    assert result["index"]["sha256"] == receipt["index"]["sha256"]
    with pytest.raises(FileExistsError):
        notices.export(*arguments, output)


@pytest.mark.parametrize("change", ["missing", "duplicate", "symlink", "traversal",
                                   "absolute", "wrong_root", "oversized", "special"])
def test_unsafe_or_incomplete_archive_rejected(tmp_path, monkeypatch, change):
    arguments = fixture(tmp_path, monkeypatch, archive_change=change)
    with pytest.raises(ValueError):
        notices.export(*arguments, tmp_path / "export")
    assert not (tmp_path / "export/SOURCE_NOTICE_INDEX.json").exists()


@pytest.mark.parametrize("change", ["receipt_pin", "archive", "source", "license",
                                   "extensions", "status", "url", "revision",
                                   "duplicate_mafft", "omitted_mafft", "wrong_notice_hash",
                                   "duplicate_fasttree", "omitted_fasttree", "source_symlink"])
def test_invalid_source_or_receipt_rejected(tmp_path, monkeypatch, change):
    arguments = fixture(tmp_path, monkeypatch)
    if change == "receipt_pin":
        arguments[1] = "0" * 64
    elif change in {"archive", "source", "license"}:
        name = {"archive": "mafft.tgz", "source": "FastTree.c", "license": "LICENSE"}[change]
        with (tmp_path / "sources" / name).open("ab") as handle:
            handle.write(b"changed")
    elif change == "source_symlink":
        path = tmp_path / "sources/FastTree.c"
        moved = tmp_path / "original"
        path.rename(moved)
        path.symlink_to(moved)
    elif change in {"extensions", "status", "url", "duplicate_mafft", "omitted_mafft", "wrong_notice_hash"}:
        def mutate(data):
            if change == "extensions": data["extensions_built"] = True
            elif change == "status": data["status"] = "failed"
            elif change == "url": data["acquisition_url"] = "https://example.invalid/other"
            elif change == "duplicate_mafft": data["source_files"].append(data["source_files"][0])
            elif change == "omitted_mafft": data["source_files"].pop()
            else: data["source_files"][0]["sha256"] = "0" * 64
        rewrite(arguments, 0, mutate)
    else:
        def mutate(data):
            if change == "revision": data["revision"] = "different"
            elif change == "duplicate_fasttree": data["acquired_files"].append(data["acquired_files"][0])
            else: data["acquired_files"].pop()
        rewrite(arguments, 2, mutate)
    with pytest.raises(ValueError):
        notices.export(*arguments, tmp_path / "export")
    assert not (tmp_path / "export/SOURCE_NOTICE_INDEX.json").exists()


@pytest.mark.parametrize("header", [b"no leading comment", b"/* unterminated", b"/*" + b"x" * 16_384 + b"*/"])
def test_invalid_header_boundaries_rejected(tmp_path, monkeypatch, header):
    arguments = fixture(tmp_path, monkeypatch, header=header)
    with pytest.raises(ValueError, match="leading notice"):
        notices.export(*arguments, tmp_path / "export")


@pytest.mark.parametrize("change", ["payload", "extra", "symlink", "directory_link",
                                   "fifo", "payload_fifo", "index_fifo"])
def test_invalid_relocated_export_rejected(tmp_path, monkeypatch, change):
    arguments = fixture(tmp_path, monkeypatch)
    output = tmp_path / "export"
    receipt = notices.export(*arguments, output)
    path = output / "fasttree/LICENSE"
    if change == "payload": path.write_bytes(b"changed")
    elif change == "extra": (output / "extra").write_bytes(b"extra")
    elif change == "symlink":
        path.unlink()
        path.symlink_to(tmp_path / "sources/LICENSE")
    elif change == "directory_link":
        moved = tmp_path / "fasttree"
        (output / "fasttree").rename(moved)
        (output / "fasttree").symlink_to(moved)
    elif change == "fifo": os.mkfifo(output / "unexpected_fifo")
    elif change == "payload_fifo":
        path.unlink()
        os.mkfifo(path)
    else:
        index = output / "SOURCE_NOTICE_INDEX.json"
        index.unlink()
        os.mkfifo(index)
    with pytest.raises(ValueError):
        notices.verify(output, receipt["index"]["sha256"])


def test_watched_input_change_prevents_finalization(tmp_path, monkeypatch):
    arguments = fixture(tmp_path, monkeypatch)
    original = notices._check
    calls = 0
    def change_late(ref):
        nonlocal calls
        calls += 1
        if calls == 11:
            with Path(ref["path"]).open("ab") as handle:
                handle.write(b"changed")
        return original(ref)
    monkeypatch.setattr(notices, "_check", change_late)
    with pytest.raises(ValueError, match="Changed external-tool"):
        notices.export(*arguments, tmp_path / "export")
    assert (tmp_path / "export").is_dir()
    assert not (tmp_path / "export/SOURCE_NOTICE_INDEX.json").exists()

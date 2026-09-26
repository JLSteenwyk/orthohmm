import hashlib
import io
import tarfile

import pytest

from benchmark_tools import build_publication_mafft as module


def archive(tmp_path, monkeypatch, members):
    path = tmp_path / "source.tgz"
    with tarfile.open(path, "w:gz") as stream:
        for name, kind in members:
            member = tarfile.TarInfo(name)
            member.type = kind
            if kind == tarfile.REGTYPE:
                member.size = 4
                stream.addfile(member, io.BytesIO(b"test"))
            else:
                member.linkname = "/tmp/forbidden"
                stream.addfile(member)
    monkeypatch.setattr(module, "SHA", hashlib.sha256(path.read_bytes()).hexdigest())
    return path


def test_regular_extract_and_overwrite_refusal(tmp_path, monkeypatch):
    source = archive(tmp_path, monkeypatch, [(module.ROOT + "/core/source.c", tarfile.REGTYPE)])
    target = tmp_path / "unpacked"
    records = module.unpack(source, target)
    assert len(records) == 1 and records[0]["bytes"] == 4
    assert (target / module.ROOT / "core/source.c").read_bytes() == b"test"
    with pytest.raises(FileExistsError):
        module.unpack(source, target)


@pytest.mark.parametrize("name,kind", [
    ("/absolute", tarfile.REGTYPE),
    (module.ROOT + "/../escape", tarfile.REGTYPE),
    ("wrong-root/file", tarfile.REGTYPE),
    (module.ROOT + "/link", tarfile.SYMTYPE),
    (module.ROOT + "/link", tarfile.LNKTYPE),
])
def test_unsafe_member_rejected(tmp_path, monkeypatch, name, kind):
    source = archive(tmp_path, monkeypatch, [(name, kind)])
    with pytest.raises(ValueError, match="Unsafe"):
        module.unpack(source, tmp_path / "unpacked")
    assert not (tmp_path / "unpacked").exists()


def test_duplicate_rejected(tmp_path, monkeypatch):
    source = archive(tmp_path, monkeypatch, [(module.ROOT + "/file", tarfile.REGTYPE)] * 2)
    with pytest.raises(ValueError, match="Duplicate"):
        module.unpack(source, tmp_path / "unpacked")


def test_wrong_archive_hash_rejected(tmp_path):
    source = tmp_path / "wrong.tgz"
    source.write_bytes(b"wrong")
    with pytest.raises(ValueError, match="SHA-256"):
        module.unpack(source, tmp_path / "unpacked")

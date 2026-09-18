import gzip
import hashlib
import io
import tarfile

import pytest

from benchmark_tools.audit_qfo_source_archive import verify_members


def archive(tmp_path, names=("selected",), symlink=False):
    path = tmp_path / "source.tar.gz"
    with tarfile.open(path, "w:gz") as stream:
        for name in names:
            member = tarfile.TarInfo(name)
            if symlink:
                member.type = tarfile.SYMTYPE
                member.linkname = "outside"
                stream.addfile(member)
            else:
                member.size = 3
                stream.addfile(member, io.BytesIO(b"ABC"))
    return path


def expected():
    return {"selected": {"bytes": 3, "sha256": hashlib.sha256(b"ABC").hexdigest()}}


def test_match_without_extraction(tmp_path):
    path = archive(tmp_path, ("other", "selected"))
    result = verify_members(path, expected())
    assert result["members"] == expected()
    assert result["tar_members_visited"] == 2
    assert result["gzip_read_to_eof"]
    assert not (tmp_path / "selected").exists()


@pytest.mark.parametrize("names,symlink", [(("other",), False), (("selected", "selected"), False), (("selected",), True)])
def test_missing_duplicate_or_symlink(tmp_path, names, symlink):
    with pytest.raises(ValueError):
        verify_members(archive(tmp_path, names, symlink), expected())


def test_wrong_digest(tmp_path):
    wanted = expected()
    wanted["selected"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="differs"):
        verify_members(archive(tmp_path), wanted)


def test_corrupt_gzip_trailer(tmp_path):
    path = archive(tmp_path)
    content = bytearray(path.read_bytes())
    content[-8] ^= 1
    path.write_bytes(content)
    with pytest.raises((gzip.BadGzipFile, tarfile.ReadError)):
        verify_members(path, expected())


def test_non_gzip_trailing_bytes_not_ignored(tmp_path):
    path = archive(tmp_path)
    path.write_bytes(path.read_bytes() + b"corrupt trailer")
    with pytest.raises(gzip.BadGzipFile):
        verify_members(path, expected())

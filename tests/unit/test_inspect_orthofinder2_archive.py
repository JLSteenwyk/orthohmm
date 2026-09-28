import hashlib
import io
import tarfile

import pytest

from benchmark_tools.inspect_orthofinder2_archive import inventory


def fixture(tmp_path):
    path = tmp_path / "example.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for name in ("README", "data/treefam2reference.txt", "data/TF.nhx", "data/other.txt"):
            member = tarfile.TarInfo(name)
            member.size = 4
            archive.addfile(member, io.BytesIO(b"test"))
        link = tarfile.TarInfo("unsafe-link")
        link.type = tarfile.SYMTYPE
        link.linkname = "/tmp/never-extracted"
        archive.addfile(link)
    return path, path.stat().st_size, hashlib.md5(path.read_bytes()).hexdigest()


def test_inventory_without_extraction(tmp_path):
    path, size, digest = fixture(tmp_path)
    result = inventory(path, size, digest)
    assert result["member_count"] == 5
    assert len(result["candidate_names"]) == 3
    assert result["originals_recovered"] is False
    assert list(tmp_path.iterdir()) == [path]


@pytest.mark.parametrize("change", ["size", "hash", "symlink"])
def test_reject_corrupt_or_indirect_archive(tmp_path, change):
    path, size, digest = fixture(tmp_path)
    if change == "size":
        size += 1
    elif change == "hash":
        digest = "0"*32
    else:
        link = tmp_path / "link"
        link.symlink_to(path)
        path = link
    with pytest.raises(ValueError):
        inventory(path, size, digest)

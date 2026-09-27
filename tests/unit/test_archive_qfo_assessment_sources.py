import hashlib
import io
import tarfile

import pytest

from benchmark_tools.archive_qfo_assessment_sources import source_inventory, verify_archive


def item(path, payload=b"x"):
    return dict(path=str(path), bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())


def test_inventory_deduplicates_equal_records(tmp_path):
    a = item(tmp_path / "a.py")
    assert source_inventory(tmp_path, [a, a, item(tmp_path / "data.txt")]) == {"sources/a.py": a}


def test_inventory_rejects_conflicting_pins(tmp_path):
    with pytest.raises(ValueError):
        source_inventory(tmp_path, [item(tmp_path / "a.py"), item(tmp_path / "a.py", b"y")])


@pytest.mark.parametrize("case", ["valid", "missing", "duplicate", "changed", "extra", "link", "traversal"])
def test_archive_readback(tmp_path, case):
    path = tmp_path / "test.tar.gz"
    names = [] if case == "missing" else ["sources/a.py"]
    if case == "duplicate":
        names *= 2
    if case == "extra":
        names.append("sources/b.py")
    if case == "traversal":
        names = ["../a.py"]
    with tarfile.open(path, "w:gz") as stream:
        for name in names:
            member = tarfile.TarInfo(name)
            if case == "link":
                member.type, member.linkname = tarfile.SYMTYPE, "/a.py"
                stream.addfile(member)
            else:
                member.size = 1
                stream.addfile(member, io.BytesIO(b"y" if case == "changed" else b"x"))
    expected = {"sources/a.py": item(tmp_path / "a.py")}
    if case == "valid":
        assert verify_archive(path, expected) == 1
    else:
        with pytest.raises(ValueError):
            verify_archive(path, expected)

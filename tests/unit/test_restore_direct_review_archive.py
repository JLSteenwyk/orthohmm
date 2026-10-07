import hashlib
import io
import json
import tarfile

import pytest

from benchmark_tools import restore_direct_review_archive as current


def fixture(tmp_path, *, index_mode=0o664, mutation=None):
    data = b"provisional review\n"
    index = dict(schema="publication_direct_review_v3", publication_ready=False,
        redistribution_clearance=False, transitive_evidence_included=False,
        files=[dict(path="review/file.txt", mode=0o644, bytes=len(data),
                    sha256=hashlib.sha256(data).hexdigest())])
    if mutation == "claim": index["publication_ready"] = True
    if mutation == "duplicate_row": index["files"].append(index["files"][0])
    raw = json.dumps(index).encode()
    entries = [(current.INDEX, raw, index_mode, tarfile.REGTYPE),
               ("review/file.txt", data, 0o644, tarfile.REGTYPE)]
    if mutation == "path": entries[1] = ("../outside", data, 0o644, tarfile.REGTYPE)
    elif mutation == "duplicate": entries.append(entries[1])
    elif mutation == "extra": entries.append(("extra", data, 0o644, tarfile.REGTYPE))
    elif mutation == "link": entries[1] = ("review/file.txt", b"", 0o644, tarfile.SYMTYPE)
    elif mutation == "mode": entries[1] = ("review/file.txt", data, 0o666, tarfile.REGTYPE)
    elif mutation == "size": entries[1] = ("review/file.txt", data + b"x", 0o644, tarfile.REGTYPE)
    elif mutation == "checksum": entries[1] = ("review/file.txt", b"X" * len(data), 0o644, tarfile.REGTYPE)
    archive = tmp_path / "review.tar.gz"
    with tarfile.open(archive, "w:gz") as stream:
        for name, content, mode, kind in entries:
            member = tarfile.TarInfo(name)
            member.size, member.mode, member.type = len(content), mode, kind
            if kind == tarfile.SYMTYPE: member.linkname = "/outside"
            stream.addfile(member, io.BytesIO(content))
    return archive, current.identity(archive)["sha256"], hashlib.sha256(raw).hexdigest(), data


@pytest.mark.parametrize("index_mode", [0o644, 0o664])
def test_actual_member_modes_and_payloads_restored_without_code_execution(tmp_path, index_mode):
    archive, sha, index_sha, data = fixture(tmp_path, index_mode=index_mode)
    output = tmp_path / "restored"
    result = current.restore(archive, sha, index_sha, output)
    assert (output / "review/file.txt").read_bytes() == data
    assert (output / current.INDEX).stat().st_mode & 0o777 == index_mode
    assert result["payloads"] == 1 and result["members"] == 2
    assert result["copied_code_executed"] is False
    with pytest.raises(FileExistsError): current.restore(archive, sha, index_sha, output)


@pytest.mark.parametrize("mutation", ["path", "duplicate", "extra", "link", "mode", "size",
                                       "checksum", "claim", "duplicate_row", "archive_anchor", "index_anchor"])
def test_invalid_archives_refused_before_any_extraction(tmp_path, mutation):
    archive, sha, index_sha, _ = fixture(tmp_path, mutation=mutation)
    if mutation == "archive_anchor": sha = "0" * 64
    elif mutation == "index_anchor": index_sha = "0" * 64
    output = tmp_path / "restored"
    with pytest.raises(ValueError): current.restore(archive, sha, index_sha, output)
    assert not output.exists()


@pytest.mark.parametrize("name", ["", ".", "./file", "/file", "a/../file", "a\\file", "a//file"])
def test_noncanonical_paths_rejected(name):
    with pytest.raises(ValueError): current.relative(name)

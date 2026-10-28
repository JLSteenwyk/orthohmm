import io
import tarfile

import pytest

from tests.archive_runtime import extract_regular_files


def archive(member):
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as handle:
        handle.addfile(member, io.BytesIO(b"data") if member.isfile() else None)
    buffer.seek(0)
    return tarfile.open(fileobj=buffer)


def test_regular_files_and_directories_are_extracted(tmp_path):
    member = tarfile.TarInfo("nested/file.txt")
    member.size = 4
    with archive(member) as handle:
        extract_regular_files(handle, tmp_path)
    assert (tmp_path / "nested/file.txt").read_bytes() == b"data"
    member = tarfile.TarInfo("empty")
    member.type = tarfile.DIRTYPE
    with archive(member) as handle:
        extract_regular_files(handle, tmp_path)
    assert (tmp_path / "empty").is_dir()


@pytest.mark.parametrize("name,kind", [("../escape", tarfile.REGTYPE),
    ("/absolute", tarfile.REGTYPE), ("symlink", tarfile.SYMTYPE),
    ("hardlink", tarfile.LNKTYPE), ("device", tarfile.CHRTYPE)])
def test_unsafe_members_are_rejected_before_extraction(tmp_path, name, kind):
    member = tarfile.TarInfo(name)
    member.type = kind
    with archive(member) as handle, pytest.raises(ValueError, match="Unsafe"):
        extract_regular_files(handle, tmp_path)
    assert not list(tmp_path.iterdir())


def test_existing_symlink_cannot_escape_destination(tmp_path):
    destination = tmp_path / "out"
    destination.mkdir()
    (destination / "nested").symlink_to(tmp_path, target_is_directory=True)
    member = tarfile.TarInfo("nested/file")
    with archive(member) as handle, pytest.raises(ValueError, match="Unsafe"):
        extract_regular_files(handle, destination)
    assert not (tmp_path / "file").exists()

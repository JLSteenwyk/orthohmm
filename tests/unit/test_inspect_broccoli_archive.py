import hashlib
import zipfile

import pytest

from benchmark_tools.inspect_broccoli_archive import inventory, run, SIZE


def fixture(tmp_path):
    path = tmp_path / "archive.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("README.txt", "fixture only")
        archive.writestr("nested/data.tar.gz", b"not inspected")
        archive.writestr("../treefam2reference.txt", "not extracted")
    return path, path.stat().st_size, hashlib.md5(path.read_bytes()).hexdigest()


def test_inventory_does_not_extract_or_claim_recovery(tmp_path):
    path, size, md5 = fixture(tmp_path)
    result = inventory(path, size, md5)
    assert result["member_count"] == 3
    assert len(result["candidate_names"]) == 2
    assert result["originals_recovered"] is False
    assert sorted(tmp_path.iterdir()) == [path]
    assert not (tmp_path.parent / "treefam2reference.txt").exists()


def test_size_mismatch(tmp_path):
    path, size, md5 = fixture(tmp_path)
    with pytest.raises(ValueError, match="size"):
        inventory(path, size + 1, md5)


def test_checksum_mismatch(tmp_path):
    path, size, _ = fixture(tmp_path)
    with pytest.raises(ValueError, match="MD5"):
        inventory(path, size, "0" * 32)


def test_symlink_rejected(tmp_path):
    path, size, md5 = fixture(tmp_path)
    link = tmp_path / "link.zip"
    link.symlink_to(path)
    with pytest.raises(ValueError, match="file type"):
        inventory(link, size, md5)


def test_resume_rejects_symlink_before_network(tmp_path):
    path, _, _ = fixture(tmp_path)
    (tmp_path / "data_Zenodo.zip").symlink_to(path)
    with pytest.raises(ValueError, match="indirect"):
        run(tmp_path)
    assert not (tmp_path / "scheduled_download_started.json").exists()


def test_resume_rejects_oversized_file_before_network(tmp_path):
    with (tmp_path / "data_Zenodo.zip").open("wb") as handle:
        handle.truncate(SIZE + 1)
    with pytest.raises(ValueError, match="exceeds"):
        run(tmp_path)
    assert not (tmp_path / "scheduled_download_started.json").exists()

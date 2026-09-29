import hashlib
from pathlib import Path
import zipfile

import pytest

from benchmark_tools.align_private_packaging import compare_payload, run


def historical():
    return [dict(path="/old/packaging/__init__.py", kind="file", bytes=3,
                 sha256=hashlib.sha256(b"old").hexdigest())]


@pytest.mark.parametrize("payload,extra", [(b"old", False), (b"new", False), (b"old", True)])
def test_exact_payload_required(tmp_path, payload, extra):
    wheel = tmp_path / "packaging.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr("packaging/__init__.py", payload)
        if extra:
            archive.writestr("packaging/extra.py", b"extra")
    if payload == b"old" and not extra:
        rows = compare_payload(wheel, historical(), Path("/old/packaging"))
        assert len(rows) == 1 and rows[0]["sha256"] == historical()[0]["sha256"]
    else:
        with pytest.raises(ValueError):
            compare_payload(wheel, historical(), Path("/old/packaging"))


def test_missing_historical_payload(tmp_path):
    with pytest.raises(ValueError, match="Missing"):
        compare_payload(tmp_path / "absent", [], Path("/old/packaging"))


def test_no_overwrite(tmp_path):
    with pytest.raises(FileExistsError):
        run(tmp_path, tmp_path / "absent", tmp_path / "absent", "0" * 64, tmp_path)


def test_changed_historical_manifest(tmp_path):
    historical_file = tmp_path / "history.json"
    historical_file.write_text("{}")
    with pytest.raises(ValueError, match="inventory changed"):
        run(tmp_path, tmp_path / "absent", historical_file, "0" * 64, tmp_path / "output")
    assert not (tmp_path / "output").exists()

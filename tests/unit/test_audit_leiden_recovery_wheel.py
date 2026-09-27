import base64
import csv
import hashlib
import io
import zipfile

import pytest

from benchmark_tools.audit_leiden_recovery_wheel import payload


def wheel(path, change=None):
    data = b"fixture"
    member = "leidenalg/__init__.py"
    if change == "escape":
        member = "../escape"
    digest = "sha256=" + base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
    name = "leidenalg-0.11.0.dist-info/RECORD"
    rows = [[member, digest, str(len(data))], [name, "", ""]]
    if change == "record_duplicate":
        rows.append(rows[0])
    if change == "bad_hash":
        rows[0][1] = "sha256=wrong"
    if change == "bad_size":
        rows[0][2] = "0"
    output = io.StringIO()
    csv.writer(output).writerows(rows)
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(member, data)
        archive.writestr(name, output.getvalue())
        if change == "extra":
            archive.writestr("leidenalg/extra.py", b"unexpected")


def test_record_verified_without_execution(tmp_path):
    path = tmp_path / "test.whl"
    wheel(path)
    assert payload(path)["leidenalg/__init__.py"] == b"fixture"


@pytest.mark.parametrize("change", ["escape", "record_duplicate", "bad_hash", "bad_size", "extra"])
def test_invalid_archives_rejected(tmp_path, change):
    path = tmp_path / "test.whl"
    wheel(path, change)
    with pytest.raises(ValueError):
        payload(path)

import hashlib
import json
import zipfile

import pytest

pytest.importorskip("sqlglot")

from benchmark_tools import inspect_selectome_tf7ab as inspector


def fixture(tmp_path, monkeypatch, members):
    archive = tmp_path / inspector.NAME
    with zipfile.ZipFile(archive, "w") as z:
        for name, data in members:
            z.writestr(name, data)
    sha = hashlib.sha256(archive.read_bytes()).hexdigest()
    monkeypatch.setattr(inspector, "SHA", sha)
    (tmp_path / "SHA256.sum").write_text(f"{sha}  {inspector.NAME}\n")
    earlier = tmp_path / "earlier.json"
    earlier.write_text(json.dumps({"inputs": [{"sha256":
        "4047be74d730e7057633b5282abfe999cc0202dd42eb194b8efc69586c2745d1"}]}))
    return archive, earlier


@pytest.mark.parametrize("members", [[], [("other.sql", "SELECT 1;")],
    [(inspector.NAME[:-4], "SELECT 1;")], [("one.sql", ""), ("two.sql", "")]])
def test_unexpected_archive_members_rejected_before_parsing(tmp_path, monkeypatch, members):
    _, earlier = fixture(tmp_path, monkeypatch, members)
    with pytest.raises(ValueError, match="archive inventory"):
        inspector.inspect(tmp_path, earlier)


def test_published_checksum_must_match(tmp_path, monkeypatch):
    _, earlier = fixture(tmp_path, monkeypatch, [])
    (tmp_path / "SHA256.sum").write_text(f"wrong  {inspector.NAME}\n")
    with pytest.raises(ValueError, match="checksum"):
        inspector.inspect(tmp_path, earlier)


def test_wrong_earlier_source_rejected(tmp_path, monkeypatch):
    _, earlier = fixture(tmp_path, monkeypatch, [])
    earlier.write_text(json.dumps({"inputs": [{"sha256": "wrong"}]}))
    with pytest.raises(ValueError, match="earlier"):
        inspector.inspect(tmp_path, earlier)

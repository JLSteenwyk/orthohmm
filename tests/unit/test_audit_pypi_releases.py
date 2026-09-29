import io
import json
from urllib.error import HTTPError

import pytest

from benchmark_tools import audit_pypi_releases as module
from benchmark_tools.audit_pypi_releases import read_pins, review


HASH = "a" * 64


def test_exact_pin():
    assert read_pins(f"# note\nBio_python==1.87 --hash=sha256:{HASH}\n")[0] == dict(
        name="bio-python", version="1.87", hashes=[HASH])


@pytest.mark.parametrize("line", ["", "pkg>=1", "pkg==1.*", "pkg[x]==1", "pkg==1",
                                 "pkg==1 --hash=md5:abc", "pkg==1 --hash=sha256:" + "z" * 64,
                                 "pkg==1;python_version<'3' --hash=sha256:" + HASH])
def test_reject_unsupported_pin(line):
    with pytest.raises(ValueError):
        read_pins(line)


def test_duplicate():
    with pytest.raises(ValueError):
        read_pins(f"pkg==1 --hash=sha256:{HASH}\npkg==2 --hash=sha256:{HASH}")


def payload():
    return dict(info=dict(name="pkg", version="1"), vulnerabilities=[
        dict(id="PYSEC-1", withdrawn=None), dict(id="PYSEC-2", withdrawn="2026-01-01")],
        urls=[dict(filename="pkg.whl", digests=dict(sha256=HASH), yanked=False)])


def test_advisories_and_artifacts():
    result = review(dict(name="pkg", version="1", hashes=[HASH]), payload())
    assert result["active_advisories"] == ["PYSEC-1"]
    assert len(result["advisories"]) == 2
    assert result["unmatched_locked_hashes"] == []


def test_custom_artifact_not_silently_matched():
    result = review(dict(name="pkg", version="1", hashes=["b" * 64]), payload())
    assert result["matching_public_artifacts"] == []
    assert result["unmatched_locked_hashes"] == ["b" * 64]


@pytest.mark.parametrize("field,value", [("info",dict(name="other",version="1")),
                                        ("info",dict(name="pkg",version="2")),
                                        ("vulnerabilities",None), ("vulnerabilities",[{}]),
                                        ("urls",None)])
def test_invalid_response(field, value):
    data = payload()
    data[field] = value
    with pytest.raises(ValueError):
        review(dict(name="pkg", version="1", hashes=[HASH]), data)


def test_missing_advisories_not_clean():
    data = payload()
    del data["vulnerabilities"]
    with pytest.raises(KeyError):
        review(dict(name="pkg", version="1", hashes=[HASH]), data)


def test_snapshot_deduplicates_queries_and_retains_response(tmp_path, monkeypatch):
    locks = [tmp_path / "first.txt", tmp_path / "second.txt"]
    for lock in locks:
        lock.write_text(f"pkg==1 --hash=sha256:{HASH}\n")
    calls = []
    def fetch(url, timeout):
        calls.append(url)
        return io.BytesIO(json.dumps(payload()).encode())
    monkeypatch.setattr(module, "urlopen", fetch)
    report = module.audit(locks, tmp_path / "out")
    assert len(calls) == report["unique_queries"] == 1
    assert len(report["rows"]) == 2
    assert report["rows"][0]["response"] == report["rows"][1]["response"]
    assert json.loads((tmp_path / "out/pkg-1.json").read_text()) == payload()
    with pytest.raises(FileExistsError):
        module.audit(locks, tmp_path / "out")


def test_http_failure_explicitly_unresolved(tmp_path, monkeypatch):
    lock = tmp_path / "lock.txt"
    lock.write_text(f"pkg==1 --hash=sha256:{HASH}\n")
    def fetch(url, timeout):
        raise HTTPError(url, 404, "Not Found", {}, None)
    monkeypatch.setattr(module, "urlopen", fetch)
    report = module.audit([lock], tmp_path / "out")
    row = report["rows"][0]
    assert row["status"] == "unresolved"
    assert "active_advisories" not in row
    assert report["comprehensive_security_clearance"] is False

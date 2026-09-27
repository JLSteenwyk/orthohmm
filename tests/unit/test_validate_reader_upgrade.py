import json

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.validate_reader_upgrade import compare_reports, verify_lock


def reports(root):
    root.mkdir()
    names = ("structure", "sequences", "events", "hierarchy")
    for name in names:
        data = dict(count=3)
        if name == "sequences":
            data["biopython_version"] = "1.86"
        if name != "structure":
            data["checked_records"] = [record(root / "structure.json")]
        (root / (name + ".json")).write_text(json.dumps(data))
    (root / "result.json").write_text(json.dumps(dict(status="verified", reports={
        n: record(root / (n + ".json")) for n in names})))


def test_report_locations_can_change(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    reports(a)
    reports(b)
    assert len(compare_reports(a, b)) == 5


@pytest.mark.parametrize("name", ["structure", "sequences", "events", "hierarchy"])
def test_scientific_changes_rejected(tmp_path, name):
    a, b = tmp_path / "a", tmp_path / "b"
    reports(a)
    reports(b)
    path = b / (name + ".json")
    data = json.loads(path.read_text())
    data["count"] = 4
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="Reader report changed"):
        compare_reports(a, b)


def test_false_report_pointer_rejected(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    reports(a)
    reports(b)
    result = json.loads((b / "result.json").read_text())
    result["reports"]["events"]["sha256"] = "wrong"
    (b / "result.json").write_text(json.dumps(result))
    with pytest.raises(ValueError, match="Changed report reference"):
        compare_reports(a, b)


def test_changed_checked_file_rejected(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    reports(a)
    reports(b)
    path = b / "events.json"
    data = json.loads(path.read_text())
    data["checked_records"][0]["sha256"] = "wrong"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        compare_reports(a, b)


def test_summary_change_rejected(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    reports(a)
    reports(b)
    path = b / "result.json"
    data = json.loads(path.read_text())
    data["status"] = "different"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="Reader report changed: result"):
        compare_reports(a, b)


def test_lock_matches():
    wheel = dict(name="BioPython", version="1.87", wheel=dict(sha256="a" * 64))
    verify_lock("# lock\nbiopython==1.87 --hash=sha256:" + "a" * 64, [wheel])


def test_version_change_is_recorded(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    reports(a)
    reports(b)
    path = b / "sequences.json"
    data = json.loads(path.read_text())
    data["biopython_version"] = "1.87"
    path.write_text(json.dumps(data))
    result = json.loads((b / "result.json").read_text())
    result["reports"]["sequences"] = record(path)
    (b / "result.json").write_text(json.dumps(result))
    assert compare_reports(a, b)[1]["biopython_versions"] == dict(previous="1.86", current="1.87")


@pytest.mark.parametrize("pin", ["biopython==1.86", "biopython>=1.87", "biopython==1.*"])
def test_wrong_or_inexact_lock_rejected(pin):
    wheel = dict(name="biopython", version="1.87", wheel=dict(sha256="a" * 64))
    with pytest.raises(ValueError):
        verify_lock(pin + " --hash=sha256:" + "a" * 64, [wheel])


def test_changed_wheel_hash_rejected():
    with pytest.raises(ValueError, match="differ from reader lock"):
        verify_lock("biopython==1.87 --hash=sha256:" + "a" * 64,
                    [dict(name="biopython", version="1.87", wheel=dict(sha256="b" * 64))])

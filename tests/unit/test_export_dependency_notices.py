import copy
import hashlib
import json

import pytest

from benchmark_tools.export_dependency_notices import safe_member, select_wheels, verify


def wheel(path="/wheel/a.whl"):
    data = b"fixture license\n"
    return dict(name="example", version="1", wheel=dict(path=path, bytes=100, sha256="a" * 64),
                notice_candidates=[dict(path="example.dist-info/LICENSE", bytes=len(data),
                                        sha256=hashlib.sha256(data).hexdigest())])


def test_shared_wheel_deduplicated_by_identity():
    inventories = [dict(status="local_wheel_notice_inventory", wheels=[wheel(p)]) for p in ("/one", "/two")]
    assert len(select_wheels(inventories)) == 1


def test_conflicting_notice_inventory_rejected():
    first, second = wheel(), wheel()
    second["notice_candidates"][0]["sha256"] = "wrong"
    with pytest.raises(ValueError, match="Conflicting inventories"):
        select_wheels([dict(status="local_wheel_notice_inventory", wheels=[first, second])])


@pytest.mark.parametrize("name", ["/escape", "../escape", "a/../b", "a//b", "a\\b"])
def test_unsafe_member_rejected(name):
    with pytest.raises(ValueError, match="Unsafe"):
        safe_member(name)


def fixture(tmp_path):
    w = wheel()
    notice = w["notice_candidates"][0]
    name = w["wheel"]["sha256"] + "/" + notice["path"]
    path = tmp_path / name
    path.parent.mkdir(parents=True)
    path.write_bytes(b"fixture license\n")
    manifest = dict(scope="wheel_notice_candidates_only", redistribution_clearance=False, wheels=[w],
        files=[dict(relative_path=name, wheel_sha256=w["wheel"]["sha256"], member=notice["path"],
                    bytes=notice["bytes"], sha256=notice["sha256"])])
    (tmp_path / "NOTICE_INDEX.json").write_text(json.dumps(manifest))
    return manifest, path


def test_export_inventory_roundtrip(tmp_path):
    fixture(tmp_path)
    assert verify(tmp_path)["files"] == 1


def test_missing_candidate_cannot_be_omitted_from_export_index(tmp_path):
    manifest, path = fixture(tmp_path)
    manifest["files"] = []
    path.unlink()
    (tmp_path / "NOTICE_INDEX.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="inventory"):
        verify(tmp_path)


def test_changed_exported_text_rejected(tmp_path):
    _, path = fixture(tmp_path)
    path.write_text("changed")
    with pytest.raises(ValueError, match="bytes differ"):
        verify(tmp_path)


def test_duplicate_wheel_inventory_rejected(tmp_path):
    manifest, _ = fixture(tmp_path)
    manifest["wheels"].append(copy.deepcopy(manifest["wheels"][0]))
    (tmp_path / "NOTICE_INDEX.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Duplicate wheel"):
        verify(tmp_path)

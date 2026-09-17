import pytest

from benchmark_tools.freeze_qfo_assessment_environment import inventory, freeze, environment, JAVA


def test_inventory_includes_nested_references(tmp_path):
    (tmp_path / "nested").mkdir()
    (tmp_path / "a").write_text("a")
    (tmp_path / "nested/b").write_text("b")
    rows = inventory(tmp_path)
    assert len(rows) == 2
    assert all(r["bytes"] == 1 and len(r["sha256"]) == 64 for r in rows)


def test_empty_inventory_rejected(tmp_path):
    with pytest.raises(ValueError):
        inventory(tmp_path)


def test_existing_snapshot_not_overwritten(tmp_path):
    with pytest.raises(FileExistsError):
        freeze(tmp_path, tmp_path)


def test_offline_java_overrides_without_parent_mutation(monkeypatch):
    monkeypatch.setenv("JAVA_HOME", "prior")
    import os
    env = environment()
    assert env["JAVA_HOME"] == str(JAVA) and env["NXF_OFFLINE"] == "true"
    assert os.environ["JAVA_HOME"] == "prior"

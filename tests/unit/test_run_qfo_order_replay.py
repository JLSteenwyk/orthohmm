import json

import pytest

from benchmark_tools import run_qfo_order_replay as module


@pytest.fixture
def plan(tmp_path):
    source = tmp_path / "source"
    source.write_text("frozen")
    path = tmp_path / "plan.json"
    value = dict(arms=list(module.ARMS), attempts_per_arm=1, profile="satellite_v2",
                 accuracy_evaluated=False, checked_records=[module.record(source)])
    module.save(path, value)
    return path, value


def test_accepts_pinned_plan(plan):
    path, value = plan
    assert module.validate_plan(path, module.record(path)["sha256"]) == value


def test_changed_plan(plan):
    path, _ = plan
    sha = module.record(path)["sha256"]
    path.write_text("{}")
    with pytest.raises(ValueError, match="Changed paired plan"):
        module.validate_plan(path, sha)


def test_changed_input(plan):
    path, value = plan
    from pathlib import Path
    Path(value["checked_records"][0]["path"]).write_text("changed")
    with pytest.raises(ValueError, match="Changed pinned record"):
        module.validate_plan(path, module.record(path)["sha256"])


@pytest.mark.parametrize("key,value", [("arms", ["canonical_order", "retained_order"]),
    ("attempts_per_arm", 2), ("profile", "satellite_v1"), ("accuracy_evaluated", True)])
def test_scope_changes_rejected(plan, key, value):
    path, data = plan
    data[key] = value
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="scope"):
        module.validate_plan(path, module.record(path)["sha256"])


def test_receipt_not_overwritten(tmp_path):
    path = tmp_path / "result.json"
    module.save(path, {"original": True})
    with pytest.raises(FileExistsError):
        module.save(path, {})

import json
from pathlib import Path

import pytest

from benchmark_tools import run_integrated_full_job as module


def plan(tmp_path):
    path = tmp_path / "plan.json"
    data = dict(dataset="orthobench", attempts=1, checkpoint_reuse=False, controlled_timing=False,
        resources=dict(cpus=32, memory_gib=128, hours=24, node="bizon"),
        expected_genes=251378, expected_refogs=70, launcher=module.record(module.__file__), checked_records=[])
    module.save(path, data)
    return path, data


def test_full_scope(tmp_path):
    path, data = plan(tmp_path)
    assert module.validate(path, module.record(path)["sha256"]) == data


@pytest.mark.parametrize("key,value", [("attempts", 2), ("checkpoint_reuse", True),
                                      ("controlled_timing", True), ("expected_genes", 16)])
def test_changed_scope(tmp_path, key, value):
    path, data = plan(tmp_path)
    data[key] = value
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="scope"):
        module.validate(path, module.record(path)["sha256"])


def test_plan_digest_required(tmp_path):
    path, _ = plan(tmp_path)
    with pytest.raises(ValueError, match="run plan"):
        module.validate(path, "wrong")


def test_changed_dependency(tmp_path):
    path, data = plan(tmp_path)
    dependency = tmp_path / "dependency"
    dependency.write_text("original")
    data["checked_records"] = [module.record(dependency)]
    path.write_text(json.dumps(data))
    dependency.write_text("changed")
    with pytest.raises(ValueError, match="frozen dependency"):
        module.validate(path, module.record(path)["sha256"])


def test_requires_allocated_job(tmp_path, monkeypatch):
    path, _ = plan(tmp_path)
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="Slurm allocation"):
        module.run(path, module.record(path)["sha256"])

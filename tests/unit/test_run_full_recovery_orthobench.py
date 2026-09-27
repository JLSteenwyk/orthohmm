import json

import pytest

from benchmark_tools import run_full_recovery_orthobench as runner


def plan(tmp_path):
    evidence = tmp_path / "evidence"
    evidence.write_text("fixed")
    return dict(attempts=1, checkpoint_reuse=False, scientific_scores_admitted=False,
                expected_genes=251378, expected_species=12,
                resource_request=dict(cpus=32, memory_gib=128, hours=24),
                source=runner.record(runner.__file__), checked_records=[runner.record(evidence)])


def write_plan(tmp_path, value):
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(value))
    return path, runner.record(path)["sha256"]


def test_valid_bound_plan(tmp_path):
    value = plan(tmp_path)
    path, sha = write_plan(tmp_path, value)
    assert runner.validate_plan(path, sha) == value


@pytest.mark.parametrize("key,value", [("attempts", 2), ("checkpoint_reuse", True),
    ("scientific_scores_admitted", True), ("expected_genes", 16), ("expected_species", 4),
    ("resource_request", dict(cpus=1, memory_gib=128, hours=24))])
def test_reject_scope_changes(tmp_path, key, value):
    data = plan(tmp_path)
    data[key] = value
    path, sha = write_plan(tmp_path, data)
    with pytest.raises(ValueError):
        runner.validate_plan(path, sha)


@pytest.mark.parametrize("fault", ["plan", "source", "evidence"])
def test_reject_changed_pins(tmp_path, fault):
    data = plan(tmp_path)
    if fault == "source":
        data["source"]["sha256"] = "wrong"
    path, sha = write_plan(tmp_path, data)
    if fault == "plan":
        sha = "wrong"
    if fault == "evidence":
        (tmp_path / "evidence").write_text("changed")
    with pytest.raises(ValueError):
        runner.validate_plan(path, sha)


def test_no_execution_without_allocation(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="Slurm"):
        runner.run(tmp_path / "absent", "unused")

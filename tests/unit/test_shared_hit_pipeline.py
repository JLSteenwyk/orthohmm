import os
from pathlib import Path

import pytest

from benchmark_tools.run_shared_hit_pipeline import selected_control
from benchmark_tools.run_membership_hillclimb import write_json


def test_control_selection_requires_complete_verified_endpoint_vector():
    row = {"dataset": "OrthoBench", "label": "canonical", "status": "independent_content_readback_verified",
           "scores": {"OrthoBench F1": .7}}
    assert selected_control({"arms": [row]}, "OrthoBench") is row
    for report in ({"arms": []}, {"arms": [row, row]}, {"arms": [{**row, "status": "running"}]}):
        with pytest.raises(ValueError, match="control unavailable"):
            selected_control(report, "OrthoBench")
    with pytest.raises(ValueError, match="endpoint"):
        selected_control({"arms": [{**row, "scores": {}}]}, "OrthoBench")


def test_worker_preserves_unverified_controls_without_running_variants(tmp_path, monkeypatch):
    import benchmark_tools.run_shared_hit_pipeline as runner
    source = Path(runner.__file__).resolve().parent.parent
    output, controls = tmp_path / "output", tmp_path / "controls"
    output.mkdir()
    controls.mkdir()
    write_json(controls / "results.json", {"source_revision": "control", "arms": []})
    plan = output / "execution_plan.json"
    write_json(plan, {"repo": str(tmp_path), "source_root": str(source), "output": str(output),
                     "source_revision": "mocked", "control_root": str(controls), "control_dependency_job": "24240",
                     "inputs": [], "sources": [], "runtime": {"files": []}, "tool_path": os.environ["PATH"],
                     "scope": "mocked control admission", "cases": [{"dataset": name, "output": str(output / name)}
                        for name in ("OrthoBench", "QfO", "ThreeKingdoms")]})
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(runner, "runtime", lambda: {"files": []})
    monkeypatch.setattr(runner, "host_activity", lambda: "mocked")
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *a, **k: "24240|FAILED|1:0\n")
    def forbidden(*args):
        raise AssertionError("Unverified controls must not launch inference")
    monkeypatch.setattr(runner, "timed_stage", forbidden)
    report = runner.run(plan)
    assert report["status"] == "finished_with_explicit_gaps"
    assert [row["status"] for row in report["arms"]] == ["not_run_unverified_control"] * 3
    assert all(not (output / name).exists() for name in ("OrthoBench", "QfO", "ThreeKingdoms"))


def test_worker_does_not_treat_a_live_control_as_terminal(tmp_path, monkeypatch):
    import benchmark_tools.run_shared_hit_pipeline as runner
    source = Path(runner.__file__).resolve().parent.parent
    output = tmp_path / "output"
    output.mkdir()
    plan = output / "execution_plan.json"
    write_json(plan, {"repo": str(tmp_path), "source_root": str(source), "output": str(output),
                     "control_dependency_job": "24240"})
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *a, **k: "24240|RUNNING|0:0\n")
    with pytest.raises(ValueError, match="authoritatively terminal"):
        runner.run(plan)
    assert not (output / "results.json").exists()

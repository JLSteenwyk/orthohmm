from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import run_singleton_tie_orthobench as runner
from benchmark_tools.run_membership_hillclimb import write_json


def control():
    return {"dataset": "OrthoBench", "label": "canonical", "status": "independent_content_readback_verified"}


def test_unique_verified_control_selected():
    row = control()
    assert runner.selected_control({"arms": [row]}) is row


@pytest.mark.parametrize("rows", [[], [control(), control()], [{**control(), "status": "failed"}],
                                 [{**control(), "label": "arrival"}]])
def test_unverified_or_ambiguous_controls_rejected(rows):
    with pytest.raises(ValueError, match="unique verified"):
        runner.selected_control({"arms": rows})


def test_failed_control_graph_never_launches_variant(tmp_path, monkeypatch):
    source = Path(runner.__file__).resolve().parent.parent
    output = tmp_path / "output"
    output.mkdir()
    path = output / "execution_plan.json"
    write_json(path, {"repo": str(tmp_path), "source_root": str(source), "output": str(output),
        "source_revision": "fixture", "inputs": [], "sources": [], "runtime": {"files": []},
        "checkpoint": str(tmp_path / "checkpoint"), "checkpoint_sha": "fixture", "case": {},
        "tool_path": "fixture", "scope": "fixture"})
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(runner, "runtime", lambda: {"files": []})
    monkeypatch.setattr(runner, "host_activity", lambda: {})
    names = [str(i) for i in range(251378)]
    species = np.arange(len(names), dtype=np.int32) % 12
    monkeypatch.setattr(runner, "load_accuracy_checkpoint", lambda *a, **k: (names, species, [], [], []))
    monkeypatch.setattr(runner, "commands", lambda *a: (["graph"], ["native"]))
    calls = []
    def failed_graph(command, *args):
        calls.append(command)
        return {"returncode": 1}
    monkeypatch.setattr(runner, "timed_stage", failed_graph)
    result = runner.run(path)
    assert len(calls) == 1 and calls[0][-2:] == ["--singleton-tie-policy", "first_hit"]
    assert [r["status"] for r in result["arms"]] == ["failed", "not_run_failed_control_gate"]
    assert not (output / "defer_ambiguous").exists()
    assert result["goal_completed"] is False and result["defaults_promoted"] is False

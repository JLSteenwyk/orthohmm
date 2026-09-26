import json
from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys

import numpy as np
import pytest

from benchmark_tools import probe_cpm_checkpoint_setup as module


def test_setup_stops_before_refinement_and_retains_mappings(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(module, "EXPECTED_COUNTS", (3, 2, 2))
    payload = tmp_path / "payload"
    payload.mkdir()
    (payload / "gene_names.txt").write_text("a\nb\nc\n")
    for name, values in (("sources", [0, 1]), ("targets", [1, 2]), ("weights", [1., 2.])):
        np.save(payload / (name + ".npy"), np.array(values))
    calls = []
    def read(path, universe):
        calls.append(path.name)
        assert universe == {"a", "b", "c"}
        return [{"a", "b"}, {"c"}]
    def forbidden(*args, **kwargs):
        pytest.fail("must not run refinement")
    replay = SimpleNamespace(
        load_replay_input=lambda **kw: (["a", "b", "c"], np.array([0, 0, 1]),
            np.array([0]), np.array([1]), np.array([2.]), {"fixture": True}),
        read_index_clusters=lambda *args: [[0, 1], [2]],
        production_refinement_hits=lambda *args: ([], [], []),
        refine_cluster_indices=forbidden)
    result = module.setup(replay, read, np, {"path": str(tmp_path / "manifest.json"), "sha256": "fixture"},
                          tmp_path, tmp_path / "final.txt")
    assert calls == ["final.txt", "orthogroups_profiles.txt", "final.txt"]
    assert (result["genes"], result["groups"], result["seed_groups"]) == (3, 2, 2)
    assert all(a["memory_mapped"] for a in result["arrays"][:3])
    assert [json.loads(line)["stage"] for line in capsys.readouterr().err.splitlines()] == [
        "before_checkpoint", "after_checkpoint", "before_checkpoint_parser", "after_checkpoint_parser",
        "after_setup", "after_setup_parser", "stopped_before_refinement"]


def prepare(tmp_path, monkeypatch):
    results = tmp_path / "benchmark_tools/results"
    results.mkdir(parents=True)
    plan = tmp_path / "plan.json"
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"files": {}}))
    plan.write_text(json.dumps({"checkpoint_manifest": module.helper.record(manifest)}))
    monkeypatch.setattr(module, "PLAN", plan.name)
    monkeypatch.setattr(module, "PLAN_SHA", module.helper.record(plan)["sha256"])
    prior = tmp_path / module.PRIOR
    prior.write_text(json.dumps(dict(checked_records=[], arms=[{}, dict(
        result={"scientific_sources": []}, environment_overrides={"PYTHONMALLOC": "debug"})])))
    monkeypatch.setattr(module, "PRIOR_SHA", module.helper.record(prior)["sha256"])
    (tmp_path / module.PROTOCOL).write_text("fixture protocol\n")
    status = tmp_path / "status.json"
    status.write_text(json.dumps({"scientific_child_command": [sys.executable]}))
    monkeypatch.setattr(module.helper, "STATUS", status.name)
    monkeypatch.setattr(module, "SETUP_PINS", {})
    return plan


@pytest.mark.parametrize("outcome", ["success", "signal", "timeout", "malformed"])
def test_one_attempt_retains_all_outcomes(tmp_path, monkeypatch, outcome):
    prepare(tmp_path, monkeypatch)
    calls = []
    def execute(command, **kwargs):
        calls.append(command)
        assert kwargs["timeout"] == 300
        assert kwargs["env"]["PYTHONMALLOC"] == "debug"
        if outcome == "timeout":
            raise subprocess.TimeoutExpired(command, 300, output=b"partial", stderr=b"boundary")
        return subprocess.CompletedProcess(command, -11 if outcome == "signal" else 0,
            stdout=b"bad" if outcome == "malformed" else b"{}", stderr=b"boundary")
    monkeypatch.setattr(module.subprocess, "run", execute)
    output = tmp_path / "out"
    if outcome == "malformed":
        with pytest.raises(json.JSONDecodeError):
            module.run(tmp_path, output)
    else:
        module.run(tmp_path, output)
    report = json.loads((output / "report.json").read_text())
    assert len(calls) == report["attempts"] == 1
    assert report["status"] == dict(success="completed", signal="failed", timeout="timed_out",
                                    malformed="checkpoint_setup_failed")[outcome]
    assert report["accuracy_admitted"] is False
    assert Path(report["stderr"]["path"]).read_bytes() == b"boundary"
    with pytest.raises(FileExistsError):
        module.run(tmp_path, output)


def test_changed_plan_refuses_launch(tmp_path, monkeypatch):
    plan = prepare(tmp_path, monkeypatch)
    plan.write_text("{}")
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: pytest.fail("unexpected launch"))
    with pytest.raises(ValueError, match="Changed import controls or checkpoint plan"):
        module.run(tmp_path, tmp_path / "out")
    assert not (tmp_path / "out").exists()

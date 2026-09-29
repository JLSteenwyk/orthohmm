import json
from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys

import numpy as np
import pytest

from benchmark_tools import probe_cpm_refinement_gc as module


def test_stage_order_and_retained_objects(tmp_path):
    payload = tmp_path / "payload"
    payload.mkdir()
    (payload / "gene_names.txt").write_text("a\nb\nc\n")
    for name, values in (("sources", [0, 1]), ("targets", [1, 2]), ("weights", [1., 2.])):
        np.save(payload / (name + ".npy"), np.array(values))
    calls = []
    clusters = [[0, 1], [2]]
    def reader(path, universe):
        calls.append("read:" + path.name)
        assert universe == {"a", "b", "c"}
        return [{"a", "b"}, {"c"}]
    def refine(received, *args, **kwargs):
        calls.append("refine")
        assert received is clusters
        assert isinstance(kwargs["rbnh_scores"], np.memmap)
        return clusters
    def write(path, received, names):
        calls.append("write")
        assert received is clusters
        path.write_text("a b\nc\n")
    replay = SimpleNamespace(
        load_replay_input=lambda **kw: (["a", "b", "c"], np.array([0, 0, 1]), [], [], [], {}),
        read_index_clusters=lambda *a: clusters, production_refinement_hits=lambda *a: ([], [], []),
        refine_cluster_indices=refine, write_clusters=write)
    result = module.stages(replay, reader, np, dict(path=str(tmp_path / "manifest"), sha256="fixture"),
                           tmp_path, tmp_path, collect_fn=lambda stage: calls.append("gc:" + stage))
    assert calls == ["read:orthogroups_profiles.txt", "gc:before_refinement", "refine",
                     "gc:after_refinement", "write", "gc:after_write", "read:refined.txt"]
    assert (result["genes"], result["species"], result["groups"]) == (3, 2, 2)


def test_collection_markers(monkeypatch, capsys):
    generations = []
    monkeypatch.setattr(module.gc, "collect", lambda generation: generations.append(generation) or 7)
    module.collect("test")
    records = [json.loads(line) for line in capsys.readouterr().err.splitlines()]
    assert generations == [2]
    assert [r["stage"] for r in records] == ["before_gc_test", "after_gc_test"]
    assert records[-1]["unreachable"] == 7


def prepare(tmp_path, monkeypatch):
    prior = tmp_path / module.PRIOR
    prior.parent.mkdir(parents=True)
    prior.write_text(json.dumps(dict(checked_records=[], result=dict(scientific_sources=[]),
        environment_overrides=dict(PYTHONMALLOC="debug"), command=[sys.executable], cwd=str(tmp_path))))
    monkeypatch.setattr(module, "PRIOR_SHA", module.helper.record(prior)["sha256"])
    protocol = tmp_path / module.PROTOCOL
    protocol.write_text("fixture\n")
    expected = tmp_path / module.setup.NATIVE / "orthogroups_profiles_refined.txt"
    expected.parent.mkdir(parents=True)
    expected.write_text("fixture\n")
    original = module.helper.record
    def record(path):
        row = original(path)
        if Path(path).name in {"orthogroups_profiles_refined.txt", "refined.txt"}:
            row["sha256"] = "f4c6f1973bc9636828baf1e6d9be3f416a18fc8ada5302fb502c082495fc1811"
        return row
    monkeypatch.setattr(module.helper, "record", record)
    monkeypatch.setattr(module.helper, "check", lambda refs: None)
    return original(protocol)["sha256"]


@pytest.mark.parametrize("outcome", ["success", "signal", "timeout", "malformed"])
def test_single_attempt_and_failure_retention(tmp_path, monkeypatch, outcome):
    protocol_sha = prepare(tmp_path, monkeypatch)
    output = tmp_path / "out"
    calls = []
    def execute(command, **kwargs):
        calls.append(command)
        assert kwargs["timeout"] == 360 and kwargs["env"]["PYTHONMALLOC"] == "debug"
        if outcome == "timeout":
            raise subprocess.TimeoutExpired(command, 360, output=b"partial", stderr=b"stage")
        (output / "refined.txt").write_text("fixture\n")
        result = dict(genes=984137, species=78, groups=390845, refinement_directed_hits=0,
            accuracy_admitted=False, affinity=[0], gc_enabled=True, gc_thresholds=[700, 10, 10],
            output=module.helper.record(output / "refined.txt"), scientific_sources=[])
        return subprocess.CompletedProcess(command, -11 if outcome == "signal" else 0,
            stdout=b"bad" if outcome == "malformed" else json.dumps(result).encode(), stderr=b"stage")
    monkeypatch.setattr(module.subprocess, "run", execute)
    if outcome == "malformed":
        with pytest.raises(json.JSONDecodeError):
            module.run(tmp_path, output, protocol_sha)
    else:
        module.run(tmp_path, output, protocol_sha)
    report = json.loads((output / "report.json").read_text())
    assert len(calls) == report["attempts"] == 1
    assert report["status"] == dict(success="completed", signal="failed", timeout="timed_out",
                                    malformed="validation_failed")[outcome]
    assert Path(report["stderr"]["path"]).read_bytes() == b"stage"
    assert report["accuracy_admitted"] is False
    with pytest.raises(FileExistsError):
        module.run(tmp_path, output, protocol_sha)


def test_changed_protocol_refuses_launch(tmp_path, monkeypatch):
    prepare(tmp_path, monkeypatch)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: pytest.fail("unexpected launch"))
    with pytest.raises(ValueError, match="protocol changed"):
        module.run(tmp_path, tmp_path / "out", "wrong")
    assert not (tmp_path / "out").exists()

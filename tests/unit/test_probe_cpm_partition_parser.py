import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import probe_cpm_partition_parser as module

ROOT = Path(__file__).resolve().parents[2]


def fixture_files(tmp_path):
    source = tmp_path / "source.py"
    source.write_text((ROOT / "benchmark_tools/audit_historical_profile_ablation.py").read_text())
    names = tmp_path / "names.txt"
    names.write_text("a\nb\nc\n")
    partition = tmp_path / "partition.txt"
    partition.write_text("a b\nc\n")
    return source, names, partition


def test_exact_function_without_module_level_execution(tmp_path):
    source, names, partition = fixture_files(tmp_path)
    source.write_text("raise RuntimeError('module must not execute')\n" + source.read_text())
    result = module.parse_only(source, names, partition)
    assert (result["genes"], result["groups"], result["memberships"]) == (3, 2, 3)
    assert result["gc_before"]["thresholds"] == result["gc_after"]["thresholds"]


@pytest.mark.parametrize("contents", ["a a\nb c\n", "a b\nb c\n", "a b\nc d\n", "a b\n"])
def test_frozen_parser_rejects_invalid_partitions(tmp_path, contents):
    source, names, partition = fixture_files(tmp_path)
    partition.write_text(contents)
    with pytest.raises(ValueError, match="partition|Partition"):
        module.parse_only(source, names, partition)


def test_duplicate_universe_and_ambiguous_function_refused(tmp_path):
    source, names, partition = fixture_files(tmp_path)
    names.write_text("a\na\n")
    with pytest.raises(ValueError, match="Duplicate"):
        module.parse_only(source, names, partition)
    source.write_text(source.read_text() + "\ndef read_partition(path, universe): pass\n")
    with pytest.raises(ValueError, match="one undecorated"):
        module.parse_only(source, names, partition)


def prepare_run(tmp_path, monkeypatch):
    paths = fixture_files(tmp_path)
    monkeypatch.setattr(module, "PINS", {p.name: module.record(p)["sha256"] for p in paths})
    status = tmp_path / "status.json"
    status.write_text(json.dumps(dict(scientific_child_command=[sys.executable],
        checked_records=[module.record(sys.executable), module.record("/usr/lib/x86_64-linux-gnu/libc.so.6")])))
    monkeypatch.setattr(module, "STATUS", status.name)
    monkeypatch.setattr(module, "STATUS_SHA", module.record(status)["sha256"])
    protocol = tmp_path / "protocol.md"
    protocol.write_text("fixture protocol\n")
    monkeypatch.setattr(module, "PROTOCOL", protocol.name)


@pytest.mark.parametrize("outcome", ["success", "signal", "timeout"])
def test_two_planned_arms_preserve_failure_without_retry(tmp_path, monkeypatch, outcome):
    prepare_run(tmp_path, monkeypatch)
    calls = []
    def execute(command, **kwargs):
        calls.append((command, kwargs))
        assert command[1:3] == ["-B", "-S"]
        assert kwargs["timeout"] == 120
        if outcome == "timeout":
            raise subprocess.TimeoutExpired(command, 120, output=b"partial", stderr=b"timeout detail")
        return subprocess.CompletedProcess(command, 0 if outcome == "success" else -11,
            stdout=b'{"fixture":true}' if outcome == "success" else b"", stderr=b"diagnostic")
    monkeypatch.setattr(module.subprocess, "run", execute)
    output = tmp_path / "out"
    report = module.run(tmp_path, output)
    assert len(calls) == 2
    assert "PYTHONMALLOC" not in calls[0][1]["env"]
    assert calls[1][1]["env"]["PYTHONMALLOC"] == "debug"
    assert report["accuracy_admitted"] is False
    for arm in report["arms"]:
        assert arm["attempts"] == 1
        assert arm["status"] == {"success": "completed", "signal": "failed", "timeout": "timed_out"}[outcome]
        assert Path(arm["stderr"]["path"]).read_bytes()
    with pytest.raises(FileExistsError):
        module.run(tmp_path, output)
    assert len(calls) == 2


def test_changed_inputs_prevent_launch(tmp_path, monkeypatch):
    prepare_run(tmp_path, monkeypatch)
    (tmp_path / "partition.txt").write_text("changed\n")
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: pytest.fail("unexpected launch"))
    with pytest.raises(ValueError, match="Changed frozen"):
        module.run(tmp_path, tmp_path / "out")
    assert not (tmp_path / "out").exists()

import json
import os
from pathlib import Path

import pytest

from orthohmm.metrics import PipelineMetrics, process_tree_rss_bytes
from orthohmm import metrics as module


@pytest.mark.skipif(not Path("/proc/self/statm").exists(), reason="Native Linux /proc RSS observation")
def test_process_tree_rss_includes_current_process():
    assert process_tree_rss_bytes(os.getpid()) > 0


def test_pipeline_metrics_writes_stage_counts_and_metadata(tmp_path):
    output = tmp_path / "metrics.json"
    with PipelineMetrics(str(output), sample_interval=0.001) as metrics:
        metrics.add_metadata(dataset="tiny")
        with metrics.stage("search"):
            payload = bytearray(1024 * 1024)
            assert payload
        metrics.add_counts(search_candidates=12, significant_hits=3)

    data = json.loads(output.read_text())
    assert data["status"] == "complete"
    assert data["metadata"]["dataset"] == "tiny"
    assert data["counts"] == {
        "search_candidates": 12,
        "significant_hits": 3,
    }
    assert data["stages"]["search"]["wall_s"] >= 0
    if Path("/proc/self/statm").exists():
        assert data["rss_measurement"] == "sampled_sum_of_linux_proc_tree_rss"
        assert data["peak_process_tree_rss_bytes"] > 0
    else:
        assert data["rss_measurement"] == "unavailable"
        assert data["peak_process_tree_rss_bytes"] == 0
        assert data["stages"]["search"]["peak_process_tree_rss_bytes"] == 0


@pytest.fixture
def missing_proc(monkeypatch):
    read_text, exists = Path.read_text, Path.exists
    def missing_read(path, *args, **kwargs):
        if str(path).startswith("/proc/"):
            raise FileNotFoundError("synthetic absent proc filesystem")
        return read_text(path, *args, **kwargs)
    def missing_exists(path):
        return False if str(path).startswith("/proc/") else exists(path)
    monkeypatch.setattr(Path, "read_text", missing_read)
    monkeypatch.setattr(Path, "exists", missing_exists)


def test_missing_proc_export_is_explicitly_unavailable(tmp_path, missing_proc):
    assert process_tree_rss_bytes(os.getpid()) == 0
    output = tmp_path / "metrics.json"
    with PipelineMetrics(str(output), sample_interval=1) as metrics:
        with metrics.stage("search"):
            metrics.add_counts(search_candidates=12)
    data = json.loads(output.read_text())
    assert data["status"] == "complete"
    assert data["rss_measurement"] == "unavailable"
    assert data["counts"] == {"search_candidates": 12}
    for row in (data, data["stages"]["search"]):
        assert row["peak_process_tree_rss_bytes"] == row["peak_process_tree_rss_gib"] == 0
        assert row["wall_s"] >= 0


def test_disabled_metrics_never_sample(monkeypatch):
    def forbidden(*args):
        pytest.fail("Disabled metrics reached resource observation")
    monkeypatch.setattr(module, "process_tree_rss_bytes", forbidden)
    monkeypatch.setattr(PipelineMetrics, "_cpu_seconds", staticmethod(forbidden))
    with PipelineMetrics(None) as metrics:
        with metrics.stage("search"):
            metrics.add_counts(genes=1)
            metrics.add_metadata(dataset="disabled")
    assert metrics.data == {} and metrics._thread is None


@pytest.mark.parametrize("code,status", [(None, "complete"), (0, "complete"), (7, "failed")])
def test_exit_status_is_retained(tmp_path, code, status):
    output = tmp_path / "metrics.json"
    with pytest.raises(SystemExit) as caught:
        with PipelineMetrics(str(output), sample_interval=1) as metrics:
            with metrics.stage("search"):
                raise SystemExit(code)
    assert caught.value.code == code
    data = json.loads(output.read_text())
    assert data["status"] == status and "search" in data["stages"]
    assert ("error" in data) == (status == "failed")
    assert metrics._thread is not None and not metrics._thread.is_alive()


def test_failed_stage_keeps_evidence_and_propagates(tmp_path):
    output = tmp_path / "metrics.json"
    with pytest.raises(RuntimeError, match="fixture failure"):
        with PipelineMetrics(str(output), sample_interval=1) as metrics:
            with metrics.stage("search"):
                metrics.add_counts(genes=1)
                raise RuntimeError("fixture failure")
    data = json.loads(output.read_text())
    assert data["status"] == "failed" and data["error"] == "RuntimeError: fixture failure"
    assert data["counts"] == {"genes": 1} and "search" in data["stages"]
    assert metrics._active_stage is None and not metrics._thread.is_alive()


@pytest.mark.parametrize("fields,expected", [("10 3", 12288), ("", 0), ("10", 0), ("10 invalid", 0)])
def test_proc_reader_converts_pages_or_preserves_invalid_read(monkeypatch, fields, expected):
    monkeypatch.setattr(Path, "read_text", lambda *_: fields)
    monkeypatch.setattr(os, "sysconf", lambda _: 4096)
    assert module._process_rss_bytes(42) == expected


@pytest.mark.parametrize("error", [FileNotFoundError, PermissionError, OSError])
def test_proc_read_failure_is_not_fabricated(monkeypatch, error):
    def fail(*args):
        raise error("synthetic proc failure")
    monkeypatch.setattr(Path, "read_text", fail)
    assert module._process_rss_bytes(42) == 0
    assert module._child_pids(42) == []


def test_process_tree_counts_each_reachable_pid_once(monkeypatch):
    rss = {42: 100, 43: 200, 44: 300}
    children = {42: [43, 44], 43: [44], 44: [42]}
    calls = []
    def read(pid):
        calls.append(pid)
        return rss[pid]
    monkeypatch.setattr(module, "_process_rss_bytes", read)
    monkeypatch.setattr(module, "_child_pids", lambda pid: children[pid])
    assert process_tree_rss_bytes(42) == 600
    assert sorted(calls) == [42, 43, 44]


def test_child_reader_ignores_non_pid_tokens(monkeypatch):
    monkeypatch.setattr(Path, "read_text", lambda _: "43 invalid 44 43")
    assert module._child_pids(42) == [43, 44, 43]

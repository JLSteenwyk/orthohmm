import copy
import json
from pathlib import Path

import pytest

from benchmark_tools.run_blast_replay_panel import validate_commands, run, query_file_records
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

ROOT = Path(__file__).resolve().parents[2]


def fixture():
    panel = json.loads((ROOT / "benchmark_tools/results/qfo_blast_replay_panel_20260923.json").read_text())
    plan = json.loads((ROOT / "benchmark_tools/results/qfo_corrected_orthomcl_prepared_20260918.json").read_text())
    return panel, plan["search_commands"]["blast"]


def test_frozen_commands():
    panel, original = fixture()
    assert validate_commands(panel, original) == Path(panel["combined"]["path"]).parent


@pytest.mark.parametrize("flag", ["-a", "-d", "-e", "-i", "-o", "-v", "-b"])
def test_changed_command_rejected(flag):
    panel, original = fixture()
    panel = copy.deepcopy(panel)
    command = panel["commands"][0]["argv"]
    command[command.index(flag) + 1] = "changed"
    with pytest.raises(ValueError, match="Changed replay"):
        validate_commands(panel, original)


def test_incomplete_panel_rejected():
    panel, original = fixture()
    panel["commands"].pop()
    with pytest.raises(ValueError):
        validate_commands(panel, original)


def test_unscheduled_execution_rejected(monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        run(Path("not_read.json"))


def test_query_metadata_does_not_invalidate_unchanged_bytes(tmp_path):
    path = tmp_path / "query.fa"
    path.write_text(">a\nAAA\n")
    item = {**record(path), "id": "a", "input_ordinal_0based": 0}
    with pytest.raises(ValueError):
        check(item)
    records = query_file_records([item])
    check(records[0])
    path.write_text(">a\nBBB\n")
    with pytest.raises(ValueError):
        check(records[0])


@pytest.mark.parametrize("change", [{"id": ""}, {"input_ordinal_0based": True},
    {"input_ordinal_0based": -1}, {"input_ordinal_0based": "0"}, {"extra": 1}])
def test_invalid_metadata_rejected(change):
    panel, _ = fixture()
    query = {**panel["queries"][0], **change}
    with pytest.raises(ValueError, match="metadata"):
        query_file_records([query])


@pytest.mark.parametrize("mutate", [False, True])
def test_runner_checks_real_query_records_before_and_after_execution(tmp_path, monkeypatch, mutate):
    from types import SimpleNamespace
    from benchmark_tools import run_blast_replay_panel as module

    panel, original = fixture()
    source = tmp_path / "source"
    source.write_text("frozen")
    for key in ("source", "plan", "runtime", "input"):
        panel[key] = record(source)
    panel["database"] = [record(source)]
    for i, query in enumerate(panel["queries"]):
        path = tmp_path / f"query_{i}.fa"
        path.write_text(f">q{i}\nAAA\n")
        query.update(record(path))
    combined = tmp_path / "combined.fa"
    combined.write_text("".join(Path(q["path"]).read_text() for q in panel["queries"]))
    panel["combined"] = record(combined)
    for command, query in zip(panel["commands"], [panel["combined"], *panel["queries"]]):
        command["argv"] = list(original)
        command["argv"][original.index("-i") + 1] = query["path"]
        command["argv"][original.index("-o") + 1] = str(tmp_path / (command["name"] + ".blast"))
    panel_path = tmp_path / "panel.json"
    panel_path.write_text(json.dumps(panel))
    monkeypatch.setenv("SLURM_JOB_ID", "test")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "180")
    monkeypatch.setattr(module, "read_frozen", lambda *a: panel)
    monkeypatch.setattr(module, "verify", lambda *a: {"search_commands": {"blast": original}, "output_root": str(tmp_path)})
    monkeypatch.setattr(module, "environment", lambda: {})
    calls = []

    def execute(argv, **kwargs):
        calls.append(argv)
        Path(argv[argv.index("-o") + 1]).write_text("diagnostic\n")
        if mutate and len(calls) == 6:
            Path(panel["queries"][0]["path"]).write_text("changed\n")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(module.subprocess, "run", execute)
    if mutate:
        with pytest.raises(ValueError, match="identity changed"):
            module.run(panel_path)
    else:
        module.run(panel_path)
    report = json.loads((tmp_path / "execution/status.json").read_text())
    assert len(calls) == 6 and len(report["stages"]) == 6
    assert report["status"] == ("failed_or_interrupted" if mutate else "native_replays_finished_pending_comparison")
    assert report["search_admitted"] is False and report["reuse_authorized"] is False
    if not mutate:
        with pytest.raises(FileExistsError):
            module.run(panel_path)

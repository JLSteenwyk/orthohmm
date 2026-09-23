import copy
import json
from pathlib import Path

import pytest

from benchmark_tools.run_blast_replay_panel import validate_commands, run

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

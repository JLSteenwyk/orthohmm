import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools.replay_lineage_read_crossing_control import replay, replay_trial

REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "benchmark_tools/results"


@pytest.fixture
def archive(tmp_path):
    target = tmp_path / "controls"
    shutil.copytree(RESULTS / "lineage_read_crossing_22018", target)
    return target


def mutate(path, change):
    data = json.loads(path.read_text())
    change(data)
    path.write_text(json.dumps(data))


@pytest.mark.parametrize("index", [0, 1, 2])
def test_real_trial_raw_replay(index, archive):
    identity = json.loads((archive / "identity.json").read_text())
    result = replay_trial(archive / f"trial_{index}", index, identity)
    assert result["result"]["control_met"] is True
    assert result["result"]["spans"]["before_to_crossing"]["root_minus_target_cpu_usec"] < 0
    assert result["result"]["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["command", "service", "removal", "nested_scope", "nested_result",
                                   "event", "outer_scope", "result"])
def test_tampered_trial_rejected(archive, fault):
    trial = archive / "trial_0"
    cases = {
        "command": ("lifecycle/command.json", lambda d: d["command"].append("unexpected")),
        "service": ("lifecycle/service.json", lambda d: d.update(returncode=1)),
        "removal": ("lifecycle/removal.json", lambda d: d["observations"][-1].update(exists=True)),
        "nested_scope": ("lifecycle/before.json", lambda d: d["target"].update(target="/other")),
        "nested_result": ("lifecycle_result.json", lambda d: d["result"].update(manager_cpu_s=0)),
        "event": ("event.json", lambda d: d.update(finished_ns=d["started_ns"])),
        "outer_scope": ("before.json", lambda d: d.update(target="/other")),
        "result": ("result.json", lambda d: d.update(control_met=False)),
    }
    filename, change = cases[fault]
    mutate(trial / filename, change)
    identity = json.loads((archive / "identity.json").read_text())
    with pytest.raises(ValueError):
        replay_trial(trial, 0, identity)


def test_complete_raw_replay_with_pinned_sources(archive):
    result = replay(archive, RESULTS / "lineage_read_crossing_scheduler_22018.txt", REPO)
    assert len(result["trials"]) == 3
    assert result["all_controls_met"] is True
    assert result["environmental_validity_established"] is False


def test_missing_trial_evidence_cannot_be_skipped(archive):
    (archive / "trial_2/result.json").unlink()
    with pytest.raises(ValueError, match="35-file"):
        replay(archive, RESULTS / "lineage_read_crossing_scheduler_22018.txt", REPO)


def test_wrong_source_hash_rejected(archive):
    mutate(archive / "identity.json", lambda d: d["sources"].update(
        {next(iter(d["sources"])): "0" * 64}))
    with pytest.raises(ValueError, match="frozen commit"):
        replay(archive, RESULTS / "lineage_read_crossing_scheduler_22018.txt", REPO)


def test_wrong_scheduler_allocation_rejected(archive, tmp_path):
    scheduler = tmp_path / "scheduler.txt"
    scheduler.write_text((RESULTS / "lineage_read_crossing_scheduler_22018.txt").read_text().replace(
        "NumCPUs=20", "NumCPUs=2"))
    with pytest.raises(ValueError, match="scheduler"):
        replay(archive, scheduler, REPO)

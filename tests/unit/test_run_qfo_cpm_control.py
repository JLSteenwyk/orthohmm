import json
from types import SimpleNamespace

import pytest

from benchmark_tools import run_qfo_cpm_control as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_run_qfo_corrected_replay import fixture, STAGES


def partitions(tmp_path, text="a b\nc\n"):
    replay, baseline = {"stages": []}, {"coverage": []}
    for i, label in enumerate(STAGES):
        old, new = tmp_path / f"old{i}.txt", tmp_path / f"new{i}.txt"
        old.write_text("a b\nc\n")
        new.write_text(text)
        replay["stages"].append({"label": label, "output": record(new)})
        baseline["coverage"].append({"label": label, "output": record(old)})
    return replay, baseline


def test_partition_equivalence_ignores_order_not_membership(tmp_path):
    replay, baseline = partitions(tmp_path, "c\nb a\n")
    rows = module.compare_partitions(replay, baseline, set("abc"))
    assert len(rows) == 4 and all(r["partition_equal"] for r in rows)
    assert not any(r["bytes_equal"] for r in rows)


def test_changed_membership_is_reported_for_each_stage(tmp_path):
    replay, baseline = partitions(tmp_path, "a\nb c\n")
    rows = module.compare_partitions(replay, baseline, set("abc"))
    assert not any(r["partition_equal"] for r in rows)


@pytest.mark.parametrize("text", ["a b\n", "a b\nb c\n", "a b\nc d\n"])
def test_invalid_control_universe_rejected(tmp_path, text):
    replay, baseline = partitions(tmp_path, text)
    with pytest.raises(ValueError):
        module.compare_partitions(replay, baseline, set("abc"))


def test_allocation_gate(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        module.run(tmp_path)


@pytest.mark.parametrize("outcome", ["equal", "different", "failed", "context", "runtime"])
def test_control_parent_preserves_outcomes_without_authorizing_variants(tmp_path, monkeypatch, outcome):
    for key, value in {"SLURM_JOB_ID": "test", "SLURM_CPUS_PER_TASK": "32", "SLURM_JOB_NODELIST": "bizon"}.items():
        monkeypatch.setenv(key, value)
    output = tmp_path / "output"
    _, baseline = partitions(tmp_path)
    plan = {"runtime": {}}
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan))
    baseline.update(status="corrected_checked_replay_admitted", plan=record(plan_path), checked_records=[])
    baseline_path = tmp_path / "benchmark_tools/results/qfo_corrected_replay_admission_21757.json"
    baseline_path.parent.mkdir(parents=True)
    baseline_path.write_text(json.dumps(baseline))
    monkeypatch.setattr(module, "BASELINE_SHA", record(baseline_path)["sha256"])
    names = tmp_path / "names.txt"
    names.write_text("a\nb\nc\n")
    context = {"output_root": str(output), "cwd": str(tmp_path), "checked_records": [],
               "environment_overrides": {"OMP_NUM_THREADS": "1"}, "expected_stages": STAGES}
    monkeypatch.setattr("benchmark_tools.checked_replay_payload_worker.corrected_evidence",
                        lambda *a: (plan, record(plan_path), {}, record(names)))
    monkeypatch.setattr("benchmark_tools.cpm_replay_context.evidence", lambda *a: context)
    checks = []
    def verify(*args):
        checks.append(args)
        return {"changed": True} if outcome == "runtime" and len(checks) > 1 else {}
    monkeypatch.setattr("benchmark_tools.verify_qfo_replay_launcher.verify", verify)
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "frozen-commit\n")
    def execute(command, **kwargs):
        assert kwargs["env"]["OMP_NUM_THREADS"] == "1"
        assert "--worker" in command
        if outcome == "failed":
            return SimpleNamespace(returncode=1)
        worker, replay = fixture()
        worker["context"] = {} if outcome == "context" else context
        for i, stage in enumerate(replay["stages"]):
            path = output / f"stage{i}.txt"
            path.write_text("a\nb c\n" if outcome == "different" else "a b\nc\n")
            stage["output"] = record(path)
        (output / "checked_worker.json").write_text(json.dumps(worker))
        (output / "replay.json").write_text(json.dumps(replay))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(module.subprocess, "run", execute)
    if outcome == "equal":
        report = module.run(tmp_path)
        assert report["status"] == "cpm_control_reproduced_pending_independent_admission"
        assert all(r["partition_equal"] for r in report["stage_comparisons"])
    else:
        with pytest.raises((ValueError, RuntimeError)):
            module.run(tmp_path)
        report = json.loads((output / "results.json").read_text())
        assert report["status"] == "failed"
        if outcome == "different":
            assert not any(r["partition_equal"] for r in report["stage_comparisons"])
    assert report["accuracy_evaluated"] is False and report["changed_arms_authorized"] is False
    with pytest.raises(FileExistsError):
        module.run(tmp_path)

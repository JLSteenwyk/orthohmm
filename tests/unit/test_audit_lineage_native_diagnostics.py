import json
from pathlib import Path

import pytest

from benchmark_tools import audit_lineage_native_diagnostics as module
from tests.unit.test_verify_lineage_native_provenance import fixture, RESULTS


@pytest.fixture
def panel(tmp_path, monkeypatch):
    archive = tmp_path / "archive"
    scheduler = tmp_path / "scheduler"
    scheduler.mkdir()
    prior = dict(runs=[])
    for index in range(3):
        context, _, prepared, verified, receipt, measured, raw = fixture(index)
        measured.update(native=dict(exit_code=0, timed_out=False, started_ns=(index * 10 + 1) * 10**9,
                                     finished_ns=(index * 10 + 3) * 10**9),
                        points=[dict(lineage=dict(boot_before="boot"))])
        verified["measurement"] = measured
        task = context["plan"]["runs"][index]
        directory = archive / Path(task["run"]["measurement_directory"]).parent.relative_to(module.ROOT)
        (directory / "measurement").mkdir(parents=True)
        for name, value in (("preparation.json", prepared), ("verification.json", verified),
                            ("lineage_diagnostic_task.json", receipt), ("measurement/lineage_report.json", measured)):
            (directory / name).write_text(json.dumps(value))
        (scheduler / f"scheduler_{module.JOBS[index]}.txt").write_text(raw)
        prior["runs"].append(dict(index=task["original_index"], method=task["method"], mode="periodic",
                                  status="validated", work_identity={"canonical": index}, evidence=[]))
    source = tmp_path / "prior.json"
    source.write_text(json.dumps(prior))
    monkeypatch.setattr(module, "prior_results", lambda path: (prior, module.record(path)))
    monkeypatch.setattr(module, "recipe_evidence", lambda *args: [])
    monkeypatch.setattr(module, "replay", lambda *args: dict(native_wall_s=2., memory={},
        original_flagged_intervals=[0], narrow_flagged_intervals=[],
        screening=dict(original_screening=dict(original_threshold_screen=dict(whole_command_screen=dict(screen_passed=True))))))
    monkeypatch.setattr(module, "validate", lambda *args: dict(input_genes=73266,
        checked_files=[], gnu_time_companion=dict(accounting={"exit_status": 0}, source=module.record(source))))
    monkeypatch.setattr(module, "fingerprint", lambda run, roots: dict(identity={"canonical": run["index"]}, evidence=[]))
    return archive, RESULTS, scheduler, source


def test_all_three_preserved_without_scientific_admission(panel):
    result = module.audit(*panel)
    assert result["validated_tasks"] == 3
    assert [r["prior_index"] for r in result["runs"]] == [1, 3, 8]
    assert all(r["original_flagged_intervals"] == [0] for r in result["runs"])
    assert result["all_outputs_equivalent"] is True
    assert result["all_narrow_intervals_pass"] is True
    assert result["temporal_issues"] == []
    assert not result["scientific_timings_admitted"]
    assert all(len(row["inventory"]) == 4 for row in result["runs"])


def test_live_job_stops_before_reading_native_archive(panel, monkeypatch):
    path = panel[2] / "scheduler_21997.txt"
    path.write_text(path.read_text().replace("JobState=COMPLETED", "JobState=RUNNING"))
    def forbidden(*args):
        pytest.fail("Read native archive before all jobs terminal")
    monkeypatch.setattr(module, "load_context", forbidden)
    with pytest.raises(ValueError, match="must be terminal"):
        module.audit(*panel)


def test_missing_scheduler_stops_before_native_read(panel, monkeypatch):
    (panel[2] / "scheduler_21997.txt").unlink()
    monkeypatch.setattr(module, "load_context", lambda *args: pytest.fail("Premature native read"))
    with pytest.raises(FileNotFoundError):
        module.audit(*panel)


def test_failed_job_is_retained(panel):
    path = panel[2] / "scheduler_21996.txt"
    path.write_text(path.read_text().replace("JobState=COMPLETED", "JobState=FAILED").replace("ExitCode=0:0", "ExitCode=1:0"))
    result = module.audit(*panel)
    assert result["validated_tasks"] == 2
    assert len(result["runs"]) == 3
    assert result["runs"][1]["status"] == "failed_or_invalid"
    assert result["all_outputs_equivalent"] is None


def test_different_output_remains_failure_with_observations(panel, monkeypatch):
    monkeypatch.setattr(module, "fingerprint", lambda *args: dict(identity={"canonical": "changed"}, evidence=[]))
    result = module.audit(*panel)
    assert result["validated_tasks"] == 0
    assert result["all_outputs_equivalent"] is False
    assert result["all_narrow_intervals_pass"] is True
    assert all(r["status"] == "output_mismatch" and r["native_wall_s"] == 2 for r in result["runs"])


def test_missing_native_report_does_not_exclude_task(panel):
    (panel[0] / "lineage_native_v1/run_01/measurement/lineage_report.json").unlink()
    result = module.audit(*panel)
    assert result["validated_tasks"] == 2
    assert result["runs"][1]["error_type"] == "FileNotFoundError"
    assert result["all_narrow_intervals_pass"] is None


def test_changed_inventory_is_detected(panel, monkeypatch):
    def change(run, roots):
        path = panel[0] / Path(run["measurement_directory"]).parent.relative_to(module.ROOT) / "unexpected"
        path.write_text("changed")
        return dict(identity={"canonical": run["index"]}, evidence=[])
    monkeypatch.setattr(module, "fingerprint", change)
    result = module.audit(*panel)
    assert all(r["status"] == "failed_or_invalid" and "inventory changed" in r["reason"] for r in result["runs"])


def test_inventory_rejects_symlinks(tmp_path):
    (tmp_path / "target").write_text("data")
    (tmp_path / "link").symlink_to(tmp_path / "target")
    with pytest.raises(ValueError, match="Symlink"):
        module.inventory(tmp_path)


def test_wrong_prior_digest_rejected(tmp_path):
    path = tmp_path / "prior.gz"
    path.write_bytes(b"not pinned")
    with pytest.raises(ValueError, match="pinned panel"):
        module.prior_results(path)


def test_replay_failure_is_not_hidden(panel, monkeypatch):
    def fail(*args):
        raise ValueError("Raw observations changed")
    monkeypatch.setattr(module, "replay", fail)
    result = module.audit(*panel)
    assert result["validated_tasks"] == 0
    assert len(result["runs"]) == 3
    assert all(r["reason"] == "Raw observations changed" for r in result["runs"])


def test_interval_flags_remain_in_report(panel, monkeypatch):
    original = module.replay
    def flagged(*args):
        result = original(*args)
        result["narrow_flagged_intervals"] = [0]
        return result
    monkeypatch.setattr(module, "replay", flagged)
    result = module.audit(*panel)
    assert result["validated_tasks"] == 3
    assert result["all_narrow_intervals_pass"] is False
    assert all(r["narrow_flagged_intervals"] == [0] for r in result["runs"])
    assert not result["environmental_validity_established"]


def test_changed_clock_domain_reported(panel):
    directory = panel[0] / "lineage_native_v1/run_01"
    path = directory / "measurement/lineage_report.json"
    measured = json.loads(path.read_text())
    measured["points"][0]["lineage"]["boot_before"] = "different"
    path.write_text(json.dumps(measured))
    path = directory / "verification.json"
    verified = json.loads(path.read_text())
    verified["measurement"] = measured
    path.write_text(json.dumps(verified))
    result = module.audit(*panel)
    assert result["temporal_issues"] == [dict(left=0, right=1), dict(left=1, right=2)]
    assert not result["scientific_timings_admitted"]

import json
from pathlib import Path

import pytest

from benchmark_tools import audit_frontier_overhead as module

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def accounting(path, failed=None, running=None, array_id=21838):
    rows = []
    for i in range(18):
        state = "FAILED" if i == failed else "RUNNING" if i == running else "COMPLETED"
        code = "1:0" if i == failed else "0:0"
        rows.append(f"{array_id}_{i}|{state}|{code}|00:10:00|20|96Gn|spark-7ff0")
    path.write_text("\n".join(rows) + "\n")


@pytest.fixture
def panel(tmp_path, monkeypatch):
    path = tmp_path / "accounting.txt"
    accounting(path)
    monkeypatch.setattr(module, "recipe_evidence", lambda *args: [])
    calls = []
    def successful(archive, context, index, scheduler):
        calls.append(index)
        task = context["plan"]["runs"][index]
        return dict(index=index, method=task["method"], mode=task["mode"], pair=task["pair"],
                    status="validated", native_wall_s=600 if task["mode"] == "boundary" else 612,
                    work_identity={"synthetic_output": task["method"]},
                    whole_command_screen_passed=True, flagged_intervals=None if task["mode"] == "boundary" else [],
                    evidence=[], boot_id="synthetic-boot", started_ns=index*1000_000_000_000,
                    finished_ns=(index*1000+612)*1_000_000_000)
    monkeypatch.setattr(module, "successful_task", successful)
    return tmp_path, path, calls, successful


def test_complete_synthetic_panel_calls_every_audit_and_never_admits(panel):
    archive, path, calls, _ = panel
    result = module.audit(archive, RESULTS, path)
    assert calls == list(range(18))
    assert result["validated_tasks"] == 18
    assert result["paired"]["numerical_budget_met"] is True
    assert result["scientific_timings_admitted"] is False
    assert result["environmental_validity_established"] is False
    assert result["publication_ready"] is False


def test_scheduler_failure_is_retained_not_audited_as_success(panel):
    archive, path, calls, _ = panel
    accounting(path, failed=4)
    result = module.audit(archive, RESULTS, path)
    assert 4 not in calls
    assert result["runs"][4]["status"] == "failed"
    assert result["failed_or_unvalidated_tasks"] == 1
    assert result["paired"]["numerical_budget_met"] is None


@pytest.mark.parametrize("error,status", [(FileNotFoundError("missing raw file"), "missing_evidence"),
                                        (ValueError("bad output"), "invalid_evidence")])
def test_audit_failure_preserved_without_hiding_other_tasks(panel, monkeypatch, error, status):
    archive, path, _, original = panel
    def fail(*args):
        if args[2] == 3:
            raise error
        return original(*args)
    monkeypatch.setattr(module, "successful_task", fail)
    result = module.audit(archive, RESULTS, path)
    assert len(result["runs"]) == 18
    assert result["runs"][3]["status"] == status
    assert result["runs"][3]["reason"] == str(error)
    assert result["validated_tasks"] == 17
    assert result["paired"]["numerical_budget_met"] is None


@pytest.mark.parametrize("panel_name,array_id", [("frontier_21838", 21838), ("pressure_21889", 21889), ("dual_21920", 21920)])
def test_nonterminal_gate_precedes_native_and_context_reads(tmp_path, monkeypatch, panel_name, array_id):
    path = tmp_path / "accounting.txt"
    accounting(path, running=17, array_id=array_id)
    def forbidden(*args):
        pytest.fail("No native/context inspection before complete terminal inventory")
    monkeypatch.setattr(module, "load_context", forbidden)
    monkeypatch.setattr(module, "successful_task", forbidden)
    with pytest.raises(ValueError, match="nonterminal"):
        module.audit(tmp_path, RESULTS, path, panel=panel_name)


@pytest.mark.parametrize("fault", ["boot", "overlap", "flags"])
def test_temporal_or_screen_flags_cannot_become_admission(panel, monkeypatch, fault):
    archive, path, _, original = panel
    def changed(*args):
        row = original(*args)
        if args[2] == 1:
            if fault == "boot":
                row["boot_id"] = "changed"
            elif fault == "overlap":
                row["started_ns"] = 0
            else:
                row["flagged_intervals"] = [3]
        return row
    monkeypatch.setattr(module, "successful_task", changed)
    result = module.audit(archive, RESULTS, path)
    assert result["paired"]["numerical_budget_met"] is True
    if fault == "flags":
        assert result["observed_screens_all_pass"] is False
    else:
        assert result["temporal_issues"]
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", [None, "changed", "missing", "extra", "symlink"])
@pytest.mark.parametrize("recipe_root", ["frontier_overhead_recipe_v1", "pressure_overhead_recipe_v2"])
def test_archived_recipe_identity(tmp_path, fault, recipe_root):
    root = tmp_path / recipe_root
    root.mkdir()
    path = root / "source.py"
    path.write_text("pass\n")
    item = module.record(path)
    recipe = {"records": [dict(item, kind="file", path=str(module.ROOT / recipe_root / "source.py"))]}
    if fault == "changed":
        path.write_text("changed\n")
    elif fault == "missing":
        path.unlink()
    elif fault == "extra":
        (root / "extra.py").write_text("pass\n")
    elif fault == "symlink":
        (root / "linked.py").symlink_to(path)
    if fault:
        with pytest.raises(ValueError):
            module.recipe_evidence(tmp_path, recipe, recipe_root)
    else:
        assert module.recipe_evidence(tmp_path, recipe, recipe_root) == [item]


@pytest.mark.parametrize("panel_name", sorted(module.PANELS))
def test_successful_task_invokes_all_independent_components(tmp_path, monkeypatch, panel_name):
    context = module.load_context(RESULTS, panel_name)
    task = context["plan"]["runs"][0]
    directory = tmp_path / Path(task["run"]["measurement_directory"]).parent.relative_to(module.ROOT)
    (directory / "measurement").mkdir(parents=True)
    measured = dict(native=dict(exit_code=0, timed_out=False, started_ns=1, finished_ns=2),
                    points=[dict(frontier=dict(boot_id="test"))])
    for name, value in (("preparation.json", {}), ("verification.json", {}), ("overhead_task.json", {}),
                        ("measurement/boundary_report.json", measured)):
        (directory / name).write_text(json.dumps(value))
    (tmp_path / "scheduler_0.txt").write_text("synthetic scheduler")
    timing = directory / "native.time.tsv"
    timing.write_text("synthetic time")
    calls = []
    def binding(*args):
        calls.append("provenance")
        return dict(job_id=123, run=task["run"], measured_argv=["/frozen"])
    def replay(*args, **kwargs):
        calls.append("replay")
        assert args[1:3] == ("boundary", 123)
        assert kwargs == ({"expected_native_pressure": True} if module.panel_spec(panel_name)["pressure"] else {})
        return dict(evidence=[], native_wall_s=1., whole_command_screen_passed=True, flagged_intervals=None,
                    memory={}, screening={"native_pressure_whole_command": {"synthetic": True}})
    def validate(*args):
        calls.append("native")
        assert args[1]["command"] == ["/frozen"]
        return dict(checked_files=[], input_genes=73266, gnu_time_companion=dict(source=module.record(timing), accounting={}))
    def fingerprint(*args):
        calls.append("canonical")
        return dict(identity={"synthetic": "same"}, evidence=[])
    for name, function in (("verify", binding), ("replay", replay), ("validate", validate), ("fingerprint", fingerprint)):
        monkeypatch.setattr(module, name, function)
    result = module.successful_task(tmp_path, context, 0, dict(allocated_cpus="20", requested_memory="96Gn", node="spark-7ff0"))
    assert calls == ["provenance", "replay", "native", "canonical"]
    assert result["status"] == "validated"
    assert result["work_identity"] == {"synthetic": "same"}
    if module.panel_spec(panel_name)["pressure"]:
        assert result["native_pressure"] == dict(whole_command={"synthetic": True}, intervals=None,
                                                diagnostic_only=True, environmental_validity_established=False)
    else:
        assert "native_pressure" not in result


@pytest.mark.parametrize("failed", [None, 4])
def test_complete_pressure_panel_retains_all_tasks_and_budget_rules(panel, failed):
    archive, path, calls, _ = panel
    accounting(path, failed=failed, array_id=21889)
    result = module.audit(archive, RESULTS, path, panel="pressure_21889")
    assert len(result["runs"]) == 18
    assert calls == [i for i in range(18) if i != failed]
    assert result["paired"]["numerical_budget_met"] is (True if failed is None else None)
    assert not result["scientific_timings_admitted"]
    assert not result["environmental_validity_established"]
    assert not result["publication_ready"]


def test_pressure_panel_rejects_legacy_accounting(panel):
    archive, path, calls, _ = panel
    with pytest.raises(ValueError):
        module.audit(archive, RESULTS, path, panel="pressure_21889")
    assert not calls


@pytest.mark.parametrize("protocol_index", [0, 1])
def test_pressure_protocol_cannot_change_after_submission(panel, protocol_index):
    archive, path, calls, _ = panel
    accounting(path, array_id=21889)
    spec = module.panel_spec("pressure_21889")
    results = archive / "results"
    results.mkdir()
    names = [spec["plan_file"], spec["recipe_file"], spec["auth_file"], *spec["protocols"]]
    for name in names:
        data = (RESULTS / name).read_bytes()
        (results / name).write_bytes(data + (b"changed" if name == spec["protocols"][protocol_index] else b""))
    with pytest.raises(ValueError, match="protocol differs"):
        module.audit(archive, results, path, panel="pressure_21889")
    assert not calls


def test_dual_adapter_preserves_original_and_narrow_flags(tmp_path, monkeypatch):
    original = dict(original_threshold_screen=dict(whole_command_screen=dict(screen_passed=True), flagged_intervals=[1, 2]))
    dual = dict(original_screening=original, narrow_flagged_intervals=[2])
    def replay(*args):
        assert args == (tmp_path, 123, ["/native"])
        return dict(screening=dual, narrow_flagged_intervals=[2], native_wall_s=600)
    monkeypatch.setattr(module, "replay_dual", replay)
    result = module.replay_task(tmp_path, {"mode": "periodic"}, 123, ["/native"], {"dual": True})
    assert result["flagged_intervals"] == [1, 2]
    assert result["narrow_flagged_intervals"] == [2]
    assert result["dual_screening"] == dual
    assert result["screening"] == original


@pytest.mark.parametrize("failed", [None, 4])
def test_complete_dual_panel_retains_failure_and_numerical_budget(panel, failed):
    from tests.unit.test_verify_frontier_overhead_provenance import fixture
    archive, path, calls, _ = panel
    accounting(path, failed=failed, array_id=21920)
    for i in range(18):
        text = fixture(i, "dual_21920")[-1]
        if i == failed:
            text = text.replace("JobState=COMPLETED", "JobState=FAILED").replace("ExitCode=0:0", "ExitCode=1:0")
        (archive / f"scheduler_{i}.txt").write_text(text)
    result = module.audit(archive, RESULTS, path, panel="dual_21920")
    assert len(result["detailed_scheduler"]) == 18
    assert len(result["runs"]) == 18
    assert calls == [i for i in range(18) if i != failed]
    assert result["paired"]["numerical_budget_met"] is (True if failed is None else None)
    assert result["scientific_timings_admitted"] is False


def test_dual_missing_detailed_record_blocks_native_reads(panel, monkeypatch):
    archive, path, calls, _ = panel
    accounting(path, array_id=21920)
    monkeypatch.setattr(module, "load_context", lambda *args: pytest.fail("Native/context read before terminal proof"))
    with pytest.raises(FileNotFoundError):
        module.audit(archive, RESULTS, path, panel="dual_21920")
    assert not calls


def test_dual_scheduler_accounting_disagreement_rejected(panel):
    from tests.unit.test_verify_frontier_overhead_provenance import fixture
    archive, path, calls, _ = panel
    accounting(path, failed=0, array_id=21920)
    (archive / "scheduler_0.txt").write_text(fixture(0, "dual_21920")[-1])
    with pytest.raises(ValueError, match="differs from accounting"):
        module.audit(archive, RESULTS, path, panel="dual_21920")
    assert not calls

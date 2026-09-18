import pytest

from benchmark_tools.run_qfo_corrected_comparator_assessment import validate_stage


def fixture():
    stage = {"status": "corrected_comparator_pairs_prepared_unscored", "accuracy_evaluated": False,
             "method": "proteinortho", "participant": "qfo_corrected_proteinortho", "job_id": "1",
             "total_pairs": 2, "retained_pairs": 2, "removed_mapping_pairs": 0,
             "pairs": {"bytes": 8, "sha256": "a"}, "filtered_pairs": {"bytes": 8, "sha256": "a"}}
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "2", "JobIDRaw": "1"}
    return stage, scheduler


def test_valid():
    stage, scheduler = fixture()
    validate_stage(stage, "proteinortho", scheduler)


@pytest.mark.parametrize("key,value", [("status", "failed"), ("accuracy_evaluated", True),
    ("method", "sonic"), ("participant", "historical"), ("job_id", "2"),
    ("total_pairs", 0), ("total_pairs", True), ("retained_pairs", 1), ("removed_mapping_pairs", 1),
    ("filtered_pairs", {"bytes": 8, "sha256": "changed"})])
def test_reject_changed_stage(key, value):
    stage, scheduler = fixture()
    stage[key] = value
    with pytest.raises(ValueError):
        validate_stage(stage, "proteinortho", scheduler)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
                                      ("NodeList", "other"), ("AllocCPUS", "1")])
def test_reject_scheduler(key, value):
    stage, scheduler = fixture()
    scheduler[key] = value
    with pytest.raises(ValueError, match="scheduler"):
        validate_stage(stage, "proteinortho", scheduler)


@pytest.mark.parametrize("exit_code", [0, 1])
def test_run_preserves_native_status(monkeypatch, tmp_path, exit_code):
    import json
    from types import SimpleNamespace
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    source = tmp_path / "source"
    source.write_text("fixture")
    output = tmp_path / "execution"
    report = {"cwd": str(output), "results": str(tmp_path / "scores"), "command": ["fixture"],
              "source": module.record(source), "verified_records": [], "environment_overrides": {},
              "accuracy_admitted": False}
    monkeypatch.setattr(module, "prepare", lambda *args: report)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=exit_code))
    if exit_code:
        with pytest.raises(RuntimeError, match="scoring failed"):
            module.run(tmp_path, "proteinortho", "fixture", 1)
    else:
        assert module.run(tmp_path, "proteinortho", "fixture", 1)["status"] == "process_succeeded_pending_independent_admission"
    final = json.loads((output / "results.json").read_text())
    assert final["exit_code"] == exit_code
    assert final["accuracy_admitted"] is False
    assert (output / "preflight.json").exists()
    assert final["status"] == ("failed" if exit_code else "process_succeeded_pending_independent_admission")


def orthofinder_fixture(method):
    from benchmark_tools.run_qfo_corrected_comparator_assessment import OF_SEMANTICS
    stage, scheduler = fixture()
    stage.update(method=method, participant="qfo_corrected_" + method,
        status="corrected_orthofinder_pairs_prepared_unscored", publication_ready=False,
        semantics=OF_SEMANTICS[method])
    return stage, scheduler


@pytest.mark.parametrize("method", ["orthofinder_full", "orthofinder_sequence_only"])
def test_orthofinder_semantics(method):
    stage, scheduler = orthofinder_fixture(method)
    validate_stage(stage, method, scheduler)


@pytest.mark.parametrize("key,value", [("status", "corrected_comparator_pairs_prepared_unscored"),
    ("semantics", "root HOG cliques"), ("publication_ready", True),
    ("method", "orthofinder_sequence_only"), ("participant", "historical_of"),
    ("retained_pairs", True), ("removed_mapping_pairs", False)])
def test_wrong_orthofinder_stage_refused(key, value):
    stage, scheduler = orthofinder_fixture("orthofinder_full")
    stage[key] = value
    with pytest.raises(ValueError):
        validate_stage(stage, "orthofinder_full", scheduler)


@pytest.mark.parametrize("method", ["orthofinder_full", "fastoma"])
def test_pending_conversion_rejected_before_reading_partial_output(tmp_path, monkeypatch, method):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n123|PENDING|0:0|0:00|bizon|2\n")
    monkeypatch.setattr(module, "read_frozen", lambda *a: pytest.fail("Partial conversion read"))
    with pytest.raises(ValueError, match="COMPLETED"):
        module.prepare(tmp_path, method, "sha", 123)


def test_all_workspaces_are_distinct():
    from benchmark_tools.run_qfo_corrected_comparator_assessment import METHODS, WORK_NAMES
    assert len(set(WORK_NAMES.values())) == len(METHODS) == 5
    assert WORK_NAMES["proteinortho"] == "qc_p"
    assert WORK_NAMES["sonic"] == "qc_s"


def test_unknown_method_rejected(tmp_path):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    with pytest.raises(ValueError, match="Unknown"):
        module.prepare(tmp_path, "unknown", "sha", 123)


def test_converter_revision_pinned(tmp_path, monkeypatch):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "wrong\n")
    with pytest.raises(ValueError, match="executor changed"):
        module.converter_source(tmp_path, "orthofinder_full")


def fastoma_fixture():
    from benchmark_tools.run_qfo_corrected_comparator_assessment import FASTOMA_SEMANTICS
    stage, scheduler = fixture()
    stage.update(status="corrected_fastoma_pairs_prepared_unscored", method="fastoma",
                 participant="qfo_corrected_fastoma", publication_ready=False,
                 semantics=FASTOMA_SEMANTICS, native_pair_rows=3, native_duplicate_relations=1)
    return stage, scheduler


def test_fastoma_stage():
    stage, scheduler = fastoma_fixture()
    validate_stage(stage, "fastoma", scheduler)


@pytest.mark.parametrize("key,value", [("status", "corrected_comparator_pairs_prepared_unscored"),
    ("semantics", "HOG cliques"), ("publication_ready", True), ("native_pair_rows", 2),
    ("native_duplicate_relations", -1), ("native_duplicate_relations", True),
    ("participant", "old_fastoma"), ("removed_mapping_pairs", 1)])
def test_bad_fastoma_stage(key, value):
    stage, scheduler = fastoma_fixture()
    stage[key] = value
    with pytest.raises(ValueError):
        validate_stage(stage, "fastoma", scheduler)


def test_fastoma_converter_revision_pinned(tmp_path, monkeypatch):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "wrong\n")
    with pytest.raises(ValueError, match="executor changed"):
        module.converter_source(tmp_path, "fastoma")


def test_fastoma_converter_source_is_frozen_worktree(tmp_path, monkeypatch):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    executor = tmp_path / "benchmarks/work/publication_qfo_corrected_fastoma_pairs_v1"
    source = executor / "benchmark_tools/prepare_qfo_corrected_fastoma_pairs.py"
    source.parent.mkdir(parents=True)
    source.write_text("fixture")
    calls = []
    monkeypatch.setattr(module.subprocess, "check_output",
                        lambda *a, **k: module.FASTOMA_CONVERTER + "\n")
    monkeypatch.setattr(module.subprocess, "run", lambda args, **k: calls.append((args, k)))
    assert module.converter_source(tmp_path, "fastoma") == module.record(source)
    assert calls == [(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--",
                      "benchmark_tools"], {"check": True})]

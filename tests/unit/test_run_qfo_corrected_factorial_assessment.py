import json
from types import SimpleNamespace

import pytest

import benchmark_tools.run_qfo_corrected_factorial_assessment as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture(index):
    kind, _, semantics = module.CONVERTERS[index % 2]
    source = {"path": "/recheck", "bytes": 1, "sha256": "check"}
    stage = {"status": f"corrected_factorial_{kind}_pairs_prepared_unscored", "cell": module.CELLS[index],
        "index": index, "participant": "ohmm_qfo_corrected_factorial_" + module.CELLS[index],
        "semantics": semantics, "accuracy_evaluated": False, "publication_ready": False, "job_id": "123",
        "total_pairs": 2, "retained_pairs": 2, "removed_mapping_pairs": 0, "expected_pairs": 2,
        "pairs": {"bytes": 8, "sha256": "pairs"}, "filtered_pairs": {"bytes": 8, "sha256": "pairs"},
        "native_admission_recheck": source, "checked_records": [source]}
    scheduler = {"JobIDRaw": "123", "State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "2"}
    return stage, scheduler


@pytest.mark.parametrize("index", range(8))
def test_all_cells(index):
    stage, scheduler = fixture(index)
    module.validate_stage(stage, index, scheduler)


@pytest.mark.parametrize("key,value", [("status", "cell_pairs_prepared_unscored"), ("cell", "wrong"),
    ("index", 0), ("participant", "ohmm_qfo_factorial_p0_c0_r1"), ("semantics", "group cliques"),
    ("accuracy_evaluated", True), ("publication_ready", True), ("job_id", "999"),
    ("total_pairs", 3), ("retained_pairs", True), ("removed_mapping_pairs", 1),
    ("filtered_pairs", {"bytes": 9, "sha256": "pairs"}), ("checked_records", [])])
def test_invalid_stage(key, value):
    stage, scheduler = fixture(1)
    stage[key] = value
    with pytest.raises(ValueError):
        module.validate_stage(stage, 1, scheduler)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
    ("AllocCPUS", "32"), ("NodeList", "spark-7ff0")])
def test_invalid_scheduler(key, value):
    stage, scheduler = fixture(0)
    scheduler[key] = value
    with pytest.raises(ValueError):
        module.validate_stage(stage, 0, scheduler)


def test_group_count_must_match():
    stage, scheduler = fixture(0)
    stage["expected_pairs"] = 3
    with pytest.raises(ValueError):
        module.validate_stage(stage, 0, scheduler)


def test_live_conversion_is_rejected_before_file_access(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n123|RUNNING|0:0|00:01:00|bizon|2\n")
    monkeypatch.setattr(module, "read_frozen", lambda *a: pytest.fail("Read partial output"))
    with pytest.raises(ValueError):
        module.prepare(tmp_path, 0, "0" * 64, "123")


@pytest.mark.parametrize("exit_code", [0, 1])
def test_execution_preserves_success_or_failure_without_admission(tmp_path, monkeypatch, exit_code):
    source = tmp_path / "source.py"
    source.write_text("# fixture\n")
    results = tmp_path / "scoring"
    results.mkdir()
    (results / "scores.json").write_text("{}")
    output = tmp_path / "execution"
    preflight = {"status": "prepared_unrun", "cwd": str(output), "results": str(results),
        "source": record(source), "verified_records": [], "environment_overrides": {"OMP_NUM_THREADS": "1"},
        "command": ["frozen-nextflow", "six-endpoints"], "accuracy_admitted": False}
    monkeypatch.setattr(module, "prepare", lambda *a: preflight.copy())
    monkeypatch.setenv("SLURM_JOB_ID", "124")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    def execute(argv, **kwargs):
        assert argv == preflight["command"]
        assert kwargs["env"]["OMP_NUM_THREADS"] == "1"
        assert kwargs["cwd"] == output
        return SimpleNamespace(returncode=exit_code)
    monkeypatch.setattr(module.subprocess, "run", execute)
    if exit_code:
        with pytest.raises(RuntimeError):
            module.run(tmp_path, 0, "hash", "123")
    else:
        module.run(tmp_path, 0, "hash", "123")
    result = json.loads((output / "results.json").read_text())
    assert result["status"] == ("failed" if exit_code else "process_succeeded_pending_independent_admission")
    assert result["exit_code"] == exit_code
    assert result["accuracy_admitted"] is False
    assert len(result["outputs"]) == 1
    assert (output / "preflight.json").exists()
    with pytest.raises(FileExistsError):
        module.run(tmp_path, 0, "hash", "123")

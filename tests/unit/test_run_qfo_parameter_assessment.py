import json
from types import SimpleNamespace

import pytest

from benchmark_tools import run_qfo_parameter_assessment as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture(index=0):
    recheck = {"path": "/recheck", "bytes": 1, "sha256": "check"}
    stage = {"status": "corrected_parameter_native_pairs_prepared_unscored", "variant": module.VARIANTS[index],
        "index": index, "participant": "ohmm_qfo_parameter_" + module.VARIANTS[index],
        "semantics": "native phylogenetically inferred pairs", "job_id": "23000", "array_task_id": str(index),
        "accuracy_evaluated": False, "publication_ready": False, "written_pairs": 2,
        "total_pairs": 2, "retained_pairs": 2, "removed_mapping_pairs": 0,
        "pairs": {"bytes": 8, "sha256": "pairs"}, "filtered_pairs": {"bytes": 8, "sha256": "pairs"},
        "native_admission_recheck": recheck, "checked_records": [recheck]}
    return stage, {"JobIDRaw": "23000"}


@pytest.mark.parametrize("index", range(4))
def test_all_four_variants(index):
    stage, scheduler = fixture(index)
    module.validate_stage(stage, index, scheduler)


@pytest.mark.parametrize("key,value", [("status", "other"), ("variant", "other"), ("index", 2),
    ("participant", "other"), ("semantics", "group cliques"), ("job_id", "999"), ("array_task_id", "1"),
    ("accuracy_evaluated", True), ("publication_ready", True), ("written_pairs", 3),
    ("total_pairs", 3), ("retained_pairs", True), ("removed_mapping_pairs", 1),
    ("filtered_pairs", {"bytes": 9, "sha256": "pairs"}), ("checked_records", [])])
def test_reject_changed_conversion(key, value):
    stage, scheduler = fixture()
    stage[key] = value
    with pytest.raises(ValueError):
        module.validate_stage(stage, 0, scheduler)


@pytest.mark.parametrize("index", [-1, 4, True, "0", 0.0])
def test_invalid_index(index):
    with pytest.raises(ValueError):
        module.validate_stage({}, index, {})
    with pytest.raises(ValueError):
        module.completed_conversion("", index)


def test_live_job_gate_precedes_artifact_read(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n22039_0|23000|RUNNING|0:0|bizon|2\n")
    monkeypatch.setattr(module, "record", lambda *a: pytest.fail("Read unfinished conversion"))
    with pytest.raises(ValueError, match="terminal"):
        module.prepare(tmp_path, 0)


def test_successful_raw_job_binding():
    row = module.completed_conversion(
        "JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n22039_0|23000|COMPLETED|0:0|bizon|2\n", 0)
    assert row["JobIDRaw"] == "23000"


def test_original_conversion_cannot_supply_replacement():
    with pytest.raises(ValueError, match="terminal"):
        module.completed_conversion(
            "JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n21939_0|23000|COMPLETED|0:0|bizon|2\n", 0)


@pytest.mark.parametrize("outcome", ["success", "exit_failure", "exception", "source_changed"])
def test_scoring_retains_evidence_without_admitting_accuracy(tmp_path, monkeypatch, outcome):
    source = tmp_path / "source.py"
    source.write_text("# frozen\n")
    results = tmp_path / "results"
    results.mkdir()
    (results / "score.json").write_text("{}")
    output = tmp_path / "execution"
    prepared = {"status": "prepared_unrun", "cwd": str(output), "results": str(results),
        "source": record(source), "verified_records": [], "environment_overrides": {"OMP_NUM_THREADS": "1"},
        "command": ["frozen-nextflow", "six-endpoints"], "accuracy_admitted": False}
    monkeypatch.setattr(module, "prepare", lambda *a: prepared.copy())
    for key, value in {"SLURM_JOB_ID": "23100", "SLURM_CPUS_PER_TASK": "8",
                       "SLURM_JOB_NODELIST": "bizon", "SLURM_ARRAY_TASK_ID": "0"}.items():
        monkeypatch.setenv(key, value)
    def execute(argv, **kwargs):
        assert argv == prepared["command"]
        assert kwargs["cwd"] == output
        assert kwargs["env"]["OMP_NUM_THREADS"] == "1"
        if outcome == "exception":
            raise RuntimeError("Interrupted")
        if outcome == "source_changed":
            source.write_text("# changed\n")
        return SimpleNamespace(returncode=1 if outcome == "exit_failure" else 0)
    monkeypatch.setattr(module.subprocess, "run", execute)
    if outcome == "success":
        module.run(tmp_path, 0)
    else:
        with pytest.raises((RuntimeError, ValueError)):
            module.run(tmp_path, 0)
    report = json.loads((output / "results.json").read_text())
    assert report["status"] == ("process_succeeded_pending_independent_admission" if outcome == "success" else "failed")
    assert report["accuracy_admitted"] is False
    assert report["array_task_id"] == "0"
    assert (output / "preflight.json").exists()
    with pytest.raises(FileExistsError):
        module.run(tmp_path, 0)


def test_allocation_gate_precedes_preparation(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    monkeypatch.setattr(module, "prepare", lambda *a: pytest.fail("Prepared without allocation"))
    with pytest.raises(ValueError, match="scheduled"):
        module.run(tmp_path, 0)

import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import run_qfo_sequence_assessment as module


def stage(count=2):
    pair = dict(path="/pairs", bytes=count * 4, sha256="pairs")
    return dict(variant="all_hits", status="corrected_sequence_group_pairs_prepared_unscored",
        participant="ohmm_qfo_corrected_sequence_all_hits", semantics="cross-species group-derived clique pairs",
        accuracy_evaluated=False, publication_ready=False, job_id="1", total_pairs=count, expected_pairs=count,
        retained_pairs=count, removed_mapping_pairs=0, pairs=pair, filtered_pairs=dict(pair, path="/filtered"),
        graph_admission={"path": "/admission"}, prediction={"path": "/groups"},
        checked_records=[{"path": "/admission"}, {"path": "/groups"}])


@pytest.mark.parametrize("problem", [None, "zero", "variant", "status", "semantics", "job", "count", "negative",
                                    "boolean", "hash", "inventory", "live"])
def test_stage_gate(problem):
    value = stage(0 if problem == "zero" else 2)
    scheduler = dict(State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2", JobIDRaw="1")
    if problem == "variant":
        value["variant"] = "top100"
    elif problem == "status":
        value["status"] = "failed"
    elif problem == "semantics":
        value["semantics"] = "native phylogenetic pairs"
    elif problem == "job":
        value["job_id"] = "2"
    elif problem == "count":
        value["expected_pairs"] = 3
    elif problem == "negative":
        value["total_pairs"] = -1
    elif problem == "boolean":
        value["removed_mapping_pairs"] = False
    elif problem == "hash":
        value["filtered_pairs"]["sha256"] = "changed"
    elif problem == "inventory":
        value["checked_records"] = []
    elif problem == "live":
        scheduler["State"] = "RUNNING"
    if problem in (None, "zero"):
        module.validate_stage(value, "all_hits", scheduler)
    else:
        with pytest.raises(ValueError):
            module.validate_stage(value, "all_hits", scheduler)


@pytest.mark.parametrize("exit_code", [0, 3])
def test_execution_success_is_not_admission(tmp_path, monkeypatch, exit_code):
    directory = tmp_path / "execution"
    report = dict(status="prepared_unrun", source=module.record(Path(module.__file__)), verified_records=[],
        cwd=str(directory), results=str(tmp_path / "native"), environment_overrides={}, accuracy_admitted=False,
        command=[sys.executable, "-c", f"raise SystemExit({exit_code})"])
    monkeypatch.setattr(module, "prepare", lambda *args: report.copy())
    monkeypatch.setenv("SLURM_JOB_ID", "fixture")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    if exit_code:
        with pytest.raises(RuntimeError):
            module.run(tmp_path, "all_hits", "hash", "1")
    else:
        module.run(tmp_path, "all_hits", "hash", "1")
    actual = json.loads((directory / "results.json").read_text())
    assert actual["status"] == ("failed" if exit_code else "process_succeeded_pending_independent_admission")
    assert actual["accuracy_admitted"] is False
    assert actual["exit_code"] == exit_code
    assert (directory / "preflight.json").exists()


def test_unscheduled_run_rejected(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        module.run(tmp_path, "all_hits", "hash", "1")

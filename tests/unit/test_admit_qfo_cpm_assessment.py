from pathlib import Path

import pytest

from benchmark_tools import admit_qfo_cpm_assessment as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def accounting(index=0, state="COMPLETED", code="0:0", node="bizon", cpus="8"):
    return ("JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n"
            f"22074_{index}|22075|{state}|{code}|{node}|{cpus}\n")


@pytest.mark.parametrize("index", range(2))
def test_completed_array_binds_raw_id(index):
    assert module.completed_assessment(accounting(index), index)["JobIDRaw"] == "22075"


@pytest.mark.parametrize("change", [{"state": "RUNNING"}, {"state": "FAILED"},
    {"code": "1:0"}, {"node": "spark-7ff0"}, {"cpus": "2"}, {"index": 1}])
def test_wrong_scheduler_rejected(change):
    with pytest.raises(ValueError):
        module.completed_assessment(accounting(**change), 0)


@pytest.mark.parametrize("index", [-1, 2, True, "0", 0.0])
def test_invalid_index(index):
    with pytest.raises(ValueError):
        module.completed_assessment(accounting(), index)
    with pytest.raises(ValueError):
        module.admit(None, index, None)


def test_live_job_rejected_before_output_access(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: accounting(state="RUNNING"))
    monkeypatch.setattr(module, "record", lambda *a: pytest.fail("Read unfinished output"))
    with pytest.raises(ValueError):
        module.admit(tmp_path, 0, tmp_path / "report.json")


@pytest.mark.parametrize("symlink", [False, True])
def test_existing_report_not_overwritten(tmp_path, monkeypatch, symlink):
    output = tmp_path / "report.json"
    if symlink:
        output.symlink_to(tmp_path / "missing")
    else:
        output.write_text("preserved")
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: pytest.fail("Read scheduler before overwrite gate"))
    with pytest.raises(FileExistsError):
        module.admit(tmp_path, 0, output)


@pytest.mark.parametrize("problem", [None, "base", "missing_cpm", "missing_context", "missing_general",
                                   "duplicate", "foreign", "changed"])
def test_exact_input_and_cpm_helper_binding(tmp_path, problem):
    helper_dir = tmp_path / "executor/benchmark_tools"
    helper_dir.mkdir(parents=True)
    helpers = []
    for name in ("run_qfo_recovered_assessment.py", "run_simulation_methods.py",
                 "prepare_qfo_corrected_comparator_pairs.py", "run_qfo_parameter_phylogeny.py",
                 "prepare_ob_candidate_neighborhood.py", "run_qfo_cpm_phylogeny.py",
                 "run_qfo_cpm_variant.py", "cpm_replay_context.py"):
        path = helper_dir / name
        path.write_text("# frozen helper\n")
        helpers.append(record(path))
    base_path = tmp_path / "input"
    base_path.write_text("frozen input\n")
    base = [record(base_path)]
    records = base + helpers
    if problem == "base":
        records[0] = {}
    elif problem in ("missing_cpm", "missing_context", "missing_general"):
        name = {"missing_cpm": "run_qfo_cpm_variant.py", "missing_context": "cpm_replay_context.py",
                "missing_general": "prepare_ob_candidate_neighborhood.py"}[problem]
        records = [r for r in records if Path(r["path"]).name != name]
    elif problem == "duplicate":
        records.append(helpers[0])
    elif problem == "foreign":
        records.append({"path": "/foreign/helper.py"})
    elif problem == "changed":
        Path(helpers[0]["path"]).write_text("# changed\n")
    if problem is None:
        assert module.bind_records({"verified_records": records}, base, helper_dir.parent) == records
    else:
        with pytest.raises(ValueError):
            module.bind_records({"verified_records": records}, base, helper_dir.parent)

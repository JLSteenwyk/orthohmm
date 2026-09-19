from pathlib import Path

import pytest

from benchmark_tools import admit_qfo_parameter_assessment as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def accounting(state="COMPLETED", code="0:0", node="bizon", cpus="8"):
    return ("JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n"
            f"21944_0|21945|{state}|{code}|{node}|{cpus}\n")


def test_completed_array_binds_raw_id():
    assert module.completed_assessment(accounting(), 0)["JobIDRaw"] == "21945"


@pytest.mark.parametrize("change", [{"state": "RUNNING"}, {"state": "FAILED"},
    {"code": "1:0"}, {"node": "spark-7ff0"}, {"cpus": "2"}])
def test_wrong_scheduler_rejected(change):
    with pytest.raises(ValueError):
        module.completed_assessment(accounting(**change), 0)


@pytest.mark.parametrize("index", [-1, 4, True, "0", 0.0])
def test_invalid_index(index):
    with pytest.raises(ValueError):
        module.completed_assessment(accounting(), index)


def test_live_job_rejected_before_output_access(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: accounting(state="RUNNING"))
    monkeypatch.setattr(module, "record", lambda *a: pytest.fail("Read unfinished output"))
    with pytest.raises(ValueError):
        module.admit(tmp_path, 0, tmp_path / "report.json")


def helper_fixture(tmp_path):
    helper_dir = tmp_path / "executor/benchmark_tools"
    helper_dir.mkdir(parents=True)
    helpers = []
    for name in ("run_qfo_recovered_assessment.py", "run_simulation_methods.py",
                 "prepare_qfo_corrected_comparator_pairs.py", "run_qfo_parameter_phylogeny.py",
                 "prepare_ob_candidate_neighborhood.py"):
        path = helper_dir / name
        path.write_text("# frozen helper\n")
        helpers.append(record(path))
    base_path = tmp_path / "input"
    base_path.write_text("frozen input\n")
    return [record(base_path)], helpers, helper_dir.parent


def test_exact_input_and_helper_binding(tmp_path):
    base, helpers, executor = helper_fixture(tmp_path)
    records = base + helpers
    assert module.bind_records({"verified_records": records}, base, executor) == records


@pytest.mark.parametrize("problem", ["base", "missing", "duplicate", "foreign", "changed"])
def test_reject_changed_helper_or_input_evidence(tmp_path, problem):
    base, helpers, executor = helper_fixture(tmp_path)
    records = base + helpers
    if problem == "base":
        records[0] = {}
    elif problem == "missing":
        records.pop()
    elif problem == "duplicate":
        records.append(helpers[0])
    elif problem == "foreign":
        records.append({"path": "/foreign/helper.py"})
    else:
        Path(helpers[0]["path"]).write_text("# changed\n")
    with pytest.raises(ValueError):
        module.bind_records({"verified_records": records}, base, executor)


@pytest.mark.parametrize("problem", [None, "added", "removed", "changed", "two_traces", "no_trace"])
def test_full_output_inventory_and_unique_trace(tmp_path, problem):
    stats = tmp_path / "stats"
    stats.mkdir()
    trace = stats / "trace_test.txt"
    trace.write_text("trace\n")
    metric = tmp_path / "metric.json"
    metric.write_text("{}")
    report = {"outputs": [record(p) for p in sorted(tmp_path.rglob("*")) if p.is_file()]}
    if problem == "added":
        (tmp_path / "extra.json").write_text("{}")
    elif problem == "removed":
        metric.unlink()
    elif problem == "changed":
        metric.write_text("changed")
    elif problem == "two_traces":
        (stats / "trace_other.txt").write_text("trace\n")
        report["outputs"] = [record(p) for p in sorted(tmp_path.rglob("*")) if p.is_file()]
    elif problem == "no_trace":
        trace.unlink()
        report["outputs"] = [record(metric)]
    if problem:
        with pytest.raises(ValueError):
            module.check_output_inventory(report, tmp_path)
    else:
        observed, found = module.check_output_inventory(report, tmp_path)
        assert observed == report["outputs"] and found == trace

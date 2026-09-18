import json
from pathlib import Path

import pytest

from benchmark_tools import admit_qfo_corrected_blast as module


def fixture(tmp_path):
    plan = {"output_root": str(tmp_path / "native"), "search_commands": {"formatdb": ["formatdb"], "blast": ["blast"]}}
    work = Path(plan["output_root"]) / "work"
    execution = Path(plan["output_root"]) / "search_execution"
    def item(path):
        return {"path": str(path), "bytes": 4, "sha256": "fixture"}
    plan_record, runtime, source = (item(tmp_path / name) for name in ("plan.json", "runtime.json", "runner.py"))
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "180",
                 "ReqMem": "900G", "JobIDRaw": "123"}
    report = {"status": "search_exited_zero_pending_query_and_database_admission", "search_admitted": False,
              "accuracy_admitted": False, "job_id": "123", "node": "bizon", "plan": plan_record,
              "runtime": runtime, "source": source, "commands": plan["search_commands"],
              "environment": module.environment(), "cwd": str(work),
              "database_files": [item(work / ("all.fa." + s)) for s in ("phr", "pin", "psq")],
              "blast_output": item(work / "all.blast"), "stages": {}}
    for index, name in enumerate(("formatdb", "blast")):
        report["stages"][name] = {"exit_code": 0, "started_epoch": index * 2 + 1,
            "finished_epoch": index * 2 + 2, "log": item(execution / (name + ".log")),
            "timing": item(execution / (name + ".time.txt"))}
    return plan, report, scheduler, plan_record, runtime, source


def test_valid_native_execution(tmp_path):
    assert len(module.validate_execution(*fixture(tmp_path))) == 8


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
    ("NodeList", "other"), ("AllocCPUS", "32"), ("ReqMem", "900Gn"), ("JobIDRaw", "124")])
def test_scheduler_binding(tmp_path, key, value):
    args = fixture(tmp_path)
    args[2][key] = value
    with pytest.raises(ValueError, match="scheduler"):
        module.validate_execution(*args)


@pytest.mark.parametrize("key,value", [("status", "running_blast"), ("search_admitted", True),
    ("accuracy_admitted", True), ("node", "other"), ("cwd", "/foreign"),
    ("source", {}), ("plan", {}), ("runtime", {}), ("commands", {}), ("environment", {})])
def test_execution_drift(tmp_path, key, value):
    args = fixture(tmp_path)
    args[1][key] = value
    with pytest.raises(ValueError):
        module.validate_execution(*args)


@pytest.mark.parametrize("key,value", [("exit_code", 1), ("exit_code", False), ("started_epoch", float("nan")),
    ("finished_epoch", float("inf")), ("started_epoch", -1), ("started_epoch", 5)])
def test_invalid_stage(tmp_path, key, value):
    args = fixture(tmp_path)
    args[1]["stages"]["formatdb"][key] = value
    with pytest.raises(ValueError):
        module.validate_execution(*args)


@pytest.mark.parametrize("problem", ["order", "missing", "log", "timing", "database", "table"])
def test_stage_order_and_artifacts(tmp_path, problem):
    args = fixture(tmp_path)
    report = args[1]
    if problem == "order":
        report["stages"]["blast"]["started_epoch"] = 1
    elif problem == "missing":
        del report["stages"]["blast"]
    elif problem in ("log", "timing"):
        report["stages"]["blast"][problem]["path"] = "/foreign"
    elif problem == "database":
        report["database_files"].pop()
    else:
        report["blast_output"]["path"] += ".partial"
    with pytest.raises(ValueError):
        module.validate_execution(*args)


def test_pending_job_is_rejected_before_file_access(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode\n123|PENDING|0:0\n")
    monkeypatch.setattr(module, "record", lambda *a: pytest.fail("Read unfinished job artifacts"))
    with pytest.raises(ValueError, match="COMPLETED"):
        module.admit(tmp_path, 123, tmp_path / "out")


@pytest.mark.parametrize("problem", [None, "database", "mutation"])
def test_orchestration_with_mocked_component_audits(tmp_path, monkeypatch, problem):
    plan, report, scheduler, _, _, _ = fixture(tmp_path)
    def write(path, text="fixture"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return module.record(path)
    results = tmp_path / "benchmark_tools/results"
    report["plan"] = write(results / "qfo_corrected_orthomcl_prepared_20260918.json")
    report["runtime"] = write(results / "qfo_corrected_legacy_blast_runtime_20260918.json")
    report["source"] = write(tmp_path / "benchmarks/work/publication_qfo_corrected_blast_v1/benchmark_tools/run_qfo_corrected_blast.py")
    for stage in report["stages"].values():
        for key in ("log", "timing"):
            stage[key] = write(Path(stage[key]["path"]))
    report["database_files"] = [write(Path(item["path"])) for item in report["database_files"]]
    report["blast_output"] = write(Path(report["blast_output"]["path"]))
    plan.update(source=report["source"], checked_records=[], prepared_inputs=[])
    status_path = Path(plan["output_root"]) / "search_execution/status.json"
    write(status_path, json.dumps(report))
    accounting = "|".join(scheduler) + "\n" + "|".join(scheduler.values()) + "\n"
    monkeypatch.setattr(module.subprocess, "check_output", lambda args, **k:
                        accounting if args[0] == "sacct" else module.EXECUTOR + "\n")
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    monkeypatch.setattr(module, "verify", lambda *a: plan)
    calls = []
    def database(fasta, runtime, output):
        calls.append("database")
        write(output / "report.json", "{}")
        return {"status": "database_exact_sequence_parity_verified" if problem != "database" else "differences",
                "content": {"input_sequences": 984137}, "outputs": [], "checked_records": []}
    def table(blast, fasta, log, output):
        calls.append("table")
        write(output, "{}")
        if problem == "mutation":
            status_path.write_text("changed")
        return {"content": {"input_proteins": 984137, "failed_queries": 7}, "checked_records": []}
    monkeypatch.setattr(module, "audit_database", database)
    monkeypatch.setattr(module, "audit_table", table)
    output = tmp_path / "admission"
    if problem:
        with pytest.raises(ValueError):
            module.admit(tmp_path, 123, output)
    else:
        admitted = module.admit(tmp_path, 123, output)
        assert admitted["search_admitted"] is True
        assert admitted["query_coverage"]["failed_queries"] == 7
    final = json.loads((output / "report.json").read_text())
    assert final["accuracy_admitted"] is False and final["downstream_execution_authorized"] is False
    assert final["search_admitted"] == (problem is None)
    assert calls == (["database"] if problem == "database" else ["database", "table"])
    with pytest.raises(FileExistsError):
        module.admit(tmp_path, 123, output)

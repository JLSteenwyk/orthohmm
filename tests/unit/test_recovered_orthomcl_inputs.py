import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import recovered_orthomcl_inputs as module


@pytest.fixture
def evidence(tmp_path):
    root = tmp_path
    base = root / "benchmarks/results/qfo_blast_recovery_bpo_v1/checkpoint"
    paths = [base / n for n in ("all.bpo", "indexes/all_bpo.idx", "indexes/all_bpo.se")]
    paths.append(root / "benchmarks/results/qfo_corrected_orthomcl_v1/work/all.gg")
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture input")
    records = [module.record(p) for p in paths]
    coverage = dict(input_proteins=984137, hsp_rows=100, distinct_directed_pairs=50, failed_queries=2)
    report = dict(status="recovered_orthomcl_bpo_checkpoint_admitted", checkpoint_admitted=True,
                  accuracy_admitted=False, publication_ready=False, downstream_execution_authorized=False,
                  native_inputs=records[:3], checked_records=records,
                  validation=dict(content=dict(input_proteins=984137, source_hsp_rows=100,
                      source_pair_blocks=50), index_validation={"queries": 5}, outputs=[]),
                  query_coverage=coverage)
    return root, report


def test_inputs(evidence):
    root, report = evidence
    result = module.admitted_inputs(report, root, report["checked_records"][-1])
    assert set(result) == {"bpo", "offsets", "ranges", "species"}
    assert report["query_coverage"]["failed_queries"] == 2


@pytest.mark.parametrize("problem", ["status", "admission", "accuracy", "publication", "authorized",
    "scope", "rows", "pairs", "missing_native", "duplicate_native", "old_path", "missing_species",
    "conflict", "empty", "native_hash"])
def test_reject_input_contract(evidence, problem):
    root, report = evidence
    report = copy.deepcopy(report)
    species_record = copy.deepcopy(report["checked_records"][-1])
    if problem == "status":
        report["status"] = "corrected_orthomcl_bpo_checkpoint_admitted"
    elif problem == "admission":
        report["checkpoint_admitted"] = False
    elif problem in {"accuracy", "publication", "authorized"}:
        key = dict(accuracy="accuracy_admitted", publication="publication_ready",
                   authorized="downstream_execution_authorized")[problem]
        report[key] = True
    elif problem == "scope":
        report["validation"]["content"]["input_proteins"] -= 1
    elif problem in {"rows", "pairs"}:
        report["validation"]["content"]["source_hsp_rows" if problem == "rows" else "source_pair_blocks"] += 1
    elif problem == "missing_native":
        report["native_inputs"].pop()
    elif problem == "duplicate_native":
        report["native_inputs"].append(report["native_inputs"][0])
    elif problem == "old_path":
        report["native_inputs"][0]["path"] = str(root / "old/all.bpo")
    elif problem == "missing_species":
        species_record["path"] = str(root / "wrong.gg")
    elif problem == "conflict":
        report["checked_records"].append({**report["checked_records"][0], "sha256": "changed"})
    elif problem == "empty":
        report["checked_records"][-1]["bytes"] = 0
    elif problem == "native_hash":
        report["native_inputs"] = copy.deepcopy(report["native_inputs"])
        report["native_inputs"][0]["sha256"] = "changed"
    with pytest.raises(ValueError):
        module.admitted_inputs(report, root, species_record)


@pytest.mark.parametrize("problem", [None, "pending", "allocation", "revision", "source", "search",
    "coverage", "changed_input", "parent", "digest", "job_id", "node", "cpus", "memory"])
def test_verify_flow(evidence, monkeypatch, problem):
    root, report = evidence
    plan_path = root / "benchmark_tools/results/qfo_corrected_orthomcl_prepared_20260918.json"
    plan_path.parent.mkdir(parents=True)
    plan_path.write_text(json.dumps(dict(prepared_inputs=[report["checked_records"][-1]])))
    monkeypatch.setattr(module, "PLAN_SHA", module.record(plan_path)["sha256"])
    report["checked_records"].pop()
    executor = root / "benchmarks/work/executor"
    source = executor / "benchmark_tools/admit_blast_recovery_bpo.py"
    source.parent.mkdir(parents=True)
    source.write_text("fixture validator")
    monkeypatch.setattr(module, "ADMITTER_SHA", module.record(source)["sha256"])
    search_path = root / "search.json"
    search_path.write_text("{}")
    search_record = module.record(search_path)
    row = dict(JobIDRaw="123", State="COMPLETED", ExitCode="0:0", Elapsed="00:00:01",
               NodeList="bizon", AllocCPUS="2", ReqMem="64G")
    report.update(source=module.record(source), scheduler=row, recovered_search=search_record)
    report.update(admission_job_id="124", node="bizon", allocated_cpus=2, memory_mib=65536)
    for case, key, value in (("job_id", "admission_job_id", "123"), ("node", "node", "other"),
                             ("cpus", "allocated_cpus", 1), ("memory", "memory_mib", 1024)):
        if problem == case:
            report[key] = value
    report["checked_records"].append(search_record)
    if problem == "source":
        report["source"] = {**report["source"], "sha256": "wrong"}
    if problem == "parent":
        report["scheduler"] = {**row, "Elapsed": "00:00:02"}
    path = root / "benchmarks/results/qfo_blast_recovery_bpo_admission_v1/report.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(report))
    digest = module.record(path)["sha256"]
    def command(argv, **kwargs):
        if argv[0] == "git":
            return "wrong" if problem == "revision" else "fixture"
        state = "RUNNING" if problem == "pending" else "COMPLETED"
        cpus = "1" if problem == "allocation" else "2"
        return f"JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n{argv[2]}|{state}|0:0|00:00:01|bizon|{cpus}|64G\n"
    monkeypatch.setattr(module.subprocess, "check_output", command)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    def search(*args):
        if problem == "search":
            raise ValueError("search not admitted")
        coverage = dict(report["query_coverage"])
        if problem == "coverage":
            coverage["failed_queries"] = 0
        if problem == "changed_input":
            Path(report["native_inputs"][0]["path"]).write_text("mutated")
        return dict(query_coverage=coverage), [], [search_record], {}, ""
    monkeypatch.setattr(module, "verify_admission", search)
    if problem:
        with pytest.raises(ValueError):
            module.verify_inputs(root, path, "wrong" if problem == "digest" else digest, 124, executor, "fixture")
    else:
        result = module.verify_inputs(root, path, digest, 124, executor, "fixture")
        assert result["status"] == "recovered_native_input_evidence_verified_no_execution"
        assert result["execution_authorized"] is False
        assert result["query_coverage"]["failed_queries"] == 2
        assert result["scheduler"]["JobIDRaw"] == "124"


def test_source_pin_and_isolated_cli(tmp_path):
    source = Path(module.__file__).with_name("admit_blast_recovery_bpo.py")
    assert module.record(source)["sha256"] == module.ADMITTER_SHA
    subprocess.run([sys.executable, "-I", "-B", module.__file__, "--help"], cwd=tmp_path,
                   check=True, capture_output=True, text=True)

import json
from pathlib import Path

import pytest

from benchmark_tools import admit_blast_recovery_bpo as module


@pytest.mark.parametrize("problem", [None, "pending", "revision", "checkpoint", "search", "coverage",
    "row_count", "source_record", "recheck", "runtime_after", "changed_input", "native_after"])
def test_recovered_admission_flow(tmp_path, monkeypatch, problem):
    root = tmp_path
    executor = root / "benchmarks/work/executor"
    source = executor / "benchmark_tools/prepare_blast_recovery_bpo.py"
    source.parent.mkdir(parents=True)
    source.write_text("fixture source")
    monkeypatch.setattr(module, "PREPARER_SHA", module.record(source)["sha256"])
    base = root / "benchmarks/results/qfo_blast_recovery_bpo_v1"
    checkpoint_dir = base / "checkpoint"
    (checkpoint_dir / "indexes").mkdir(parents=True)
    native_inputs = []
    for name in ("all.bpo", "indexes/all_bpo.idx", "indexes/all_bpo.se"):
        path = checkpoint_dir / name
        path.write_text("fixture native data")
        native_inputs.append(module.record(path))
    search_path = root / "search.json"
    search_path.write_text("{}")
    search_record = module.record(search_path)
    content = dict(input_proteins=984137, source_hsp_rows=12, source_pair_blocks=6)
    coverage = dict(hsp_rows=12, distinct_directed_pairs=6, failed_queries=2)
    admission_scheduler = dict(JobID="22151", State="COMPLETED")
    checkpoint = dict(content=content, index_validation={"queries": 3}, checked_records=[], outputs=native_inputs)
    checkpoint_path = checkpoint_dir / "report.json"
    checkpoint_path.write_text(json.dumps(checkpoint))
    runtime = dict(status="dedicated_bpo_python_runtime_verified", manifest={"runtime": "fixture"}, mapped_files=[])
    parent = dict(status="recovered_bpo_prepared_pending_admission", source=module.record(source),
        checkpoint=module.record(checkpoint_path), job_id="123", started_epoch=1, finished_epoch=2,
        runtime_before=runtime, runtime_after=runtime, recovered_search=search_record, admission_scheduler=admission_scheduler,
        query_coverage=coverage, content=content, index_validation=checkpoint["index_validation"],
        checked_records=[search_record], accuracy_admitted=False, publication_ready=False, downstream_execution_authorized=False)
    if problem == "coverage":
        parent["query_coverage"] = {**coverage, "failed_queries": 0}
    elif problem == "row_count":
        coverage["hsp_rows"] = 99
    elif problem == "source_record":
        parent["checked_records"] = []
    (base / "report.json").write_text(json.dumps(parent))
    def output(command, **kwargs):
        if command[0] == "git":
            return "wrong" if problem == "revision" else "fixture"
        state = "PENDING" if problem == "pending" else "COMPLETED"
        return f"JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n123|{state}|0:0|00:00:01|bizon|2|64G\n"
    monkeypatch.setattr(module.subprocess, "check_output", output)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    calls = []
    def runtime_check(*args):
        calls.append("runtime")
        if problem == "runtime_after" and calls.count("runtime") == 2:
            raise ValueError("runtime changed")
        return runtime
    monkeypatch.setattr(module, "verify_runtime", runtime_check)
    def checkpoint_check(*args):
        if problem == "checkpoint":
            raise ValueError("bad checkpoint")
    monkeypatch.setattr(module, "validate_checkpoint", checkpoint_check)
    def search(*args):
        if problem == "search":
            raise ValueError("bad search")
        return dict(query_coverage=coverage), [search_record, search_record], [search_record], admission_scheduler, "fixture"
    monkeypatch.setattr(module, "verify_admission", search)
    original_read = module.read_frozen
    def read(path, digest):
        if path.name in {"qfo_corrected_orthomcl_perl_runtime_20260918.json", "qfo_corrected_orthomcl_system_helpers_20260918.json"}:
            return {}
        return original_read(path, digest)
    monkeypatch.setattr(module, "read_frozen", read)
    def native_check(*args):
        calls.append("native")
        if problem == "native_after" and calls.count("native") > 2:
            raise ValueError("native runtime changed")
    monkeypatch.setattr(module, "verify", native_check)
    def recheck(blast, fasta, directory, destination, expected_content, indexes):
        calls.append("recheck")
        assert expected_content == content
        if problem == "recheck":
            raise ValueError("content disagreement")
        destination.mkdir()
        proof = destination / "proof.json"
        proof.write_text("{}")
        if problem == "changed_input":
            search_path.write_text("mutated")
        return dict(content=content, index_validation=indexes, outputs=[module.record(proof)])
    monkeypatch.setattr(module, "recheck_content", recheck)
    destination = root / "admission"
    if problem:
        with pytest.raises(ValueError):
            module.admit(root, 123, executor, "fixture", destination)
        if problem in {"recheck", "runtime_after", "changed_input", "native_after"}:
            result = json.loads((destination / "report.json").read_text())
            assert result["status"] == "recovered_bpo_validation_failed"
            assert result["checkpoint_admitted"] is False
            assert result["downstream_execution_authorized"] is False
        else:
            assert not destination.exists()
            assert "recheck" not in calls
    else:
        result = module.admit(root, 123, executor, "fixture", destination)
        assert result == json.loads((destination / "report.json").read_text())
        assert result["checkpoint_admitted"] is True
        assert result["accuracy_admitted"] is False
        assert result["downstream_execution_authorized"] is False
        assert result["native_inputs"] == native_inputs
        assert result["query_coverage"]["failed_queries"] == 2
        assert calls.count("recheck") == 1
        with pytest.raises(FileExistsError):
            module.admit(root, 123, executor, "fixture", destination)

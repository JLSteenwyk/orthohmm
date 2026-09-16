import json
from pathlib import Path

import pytest

from benchmark_tools.assemble_simulation_results import admit_method, dataset_records, input_universe, terminal_tasks, verify_execution_status, merge_runtime_rows
from benchmark_tools.benchmark_production import file_record
from benchmark_tools.summarize_simulation_panel import METHODS
from benchmark_tools.validate_simulation_outputs import NativeOutputFailure


def test_all_tasks_must_be_uniquely_terminal():
    header = "JobID|JobIDRaw|State|ExitCode|Elapsed\n"
    first = "10_0|11|COMPLETED|0:0|00:01\n"
    second = "10_1|12|FAILED|1:0|00:01\n"
    assert len(terminal_tasks(header + first + second, 10, 2)) == 2
    for invalid in (header + first, header + first + second + second,
                    header + first + second.replace("FAILED", "RUNNING")):
        with pytest.raises(ValueError, match="not uniquely terminal"):
            terminal_tasks(invalid, 10, 2)


def test_native_failure_is_not_a_score_and_integrity_error_is_not_a_failure(monkeypatch):
    config = {"argv": ["tool"]}
    dataset = {"methods": {"orthofinder_full": config}}
    status = {"methods": {"orthofinder_full": {"status": "process_succeeded", "exit_code": 0, "argv": ["tool"]}}}
    def native_failure(*args):
        raise NativeOutputFailure("native stopped")
    monkeypatch.setattr("benchmark_tools.assemble_simulation_results.validate_orthofinder", native_failure)
    for method in ("orthofinder_full", "orthofinder_sequence_only"):
        row = admit_method(method, dataset, status, {"inputs": []}, {})
        assert row["status"] == "failed" and row["failure_stage"] == "native_output"
        assert "score" not in row
    def integrity_error(*args):
        raise ValueError("artifact changed")
    monkeypatch.setattr("benchmark_tools.assemble_simulation_results.validate_orthofinder", integrity_error)
    with pytest.raises(ValueError, match="artifact changed"):
        admit_method("orthofinder_full", dataset, status, {"inputs": []}, {})


def test_execution_and_interruption_retain_distinct_reasons():
    dataset = {"methods": {"orthofinder_full": {"argv": ["tool"]}}}
    status = {"methods": {}}
    assert admit_method("orthofinder_full", dataset, status, {}, {})["failure_stage"] == "execution_interrupted"
    status["methods"]["orthofinder_full"] = {"status": "failed", "argv": ["tool"], "exit_code": 7}
    row = admit_method("orthofinder_full", dataset, status, {}, {})
    assert row["exit_code"] == 7 and row["failure_stage"] == "execution"


def test_source_only_hmm_is_not_admitted():
    name = "orthohmm_high_sensitivity"
    dataset = {"methods": {name: {"argv": ["tool"]}}}
    status = {"methods": {name: {"status": "process_succeeded", "argv": ["tool"], "exit_code": 0}}}
    row = admit_method(name, dataset, status, {}, {})
    assert row["failure_stage"] == "native_runtime_unverified"
    assert "score" not in row
    with pytest.raises(ValueError, match="per-method native runtime"):
        admit_method(name, dataset, status, {}, {"native_runtime": {"sha256": "expected"}})


def test_merge_preserves_original_comparator_failure_and_provenance():
    rows = [{"method": m, "condition": "baseline", "seed": 1, "truth_sha256": "truth",
             "status": "complete", "scheduler": "new"} for m in METHODS]
    old = [{**r, "status": "failed", "scheduler": "old"} for r in rows]
    result = merge_runtime_rows(rows, old)
    for row in result:
        reused = row["method"].startswith("orthofinder_")
        assert row["reused_comparator"] == reused
        assert row["scheduler"] == ("old" if reused else "new")
        assert row["status"] == ("failed" if reused else "complete")
    old[0]["truth_sha256"] = "different"
    with pytest.raises(ValueError, match="different datasets or truth"):
        merge_runtime_rows(rows, old)


def test_input_universe_checks_duplicates_and_truth_counts(tmp_path):
    path = tmp_path / "species.fasta"
    path.write_text(">a\nAAA\n>a\nCCC\n")
    records = [dict(file_record(path, tmp_path), absolute_path=str(path))]
    truth = {"extant_genes": 2, "species": ["species"], "ortholog_pairs": [], "ortholog_pair_count": 0}
    with pytest.raises(ValueError, match="Duplicate"):
        input_universe(records, truth)


def test_execution_status_is_bound_to_task_inputs_and_sources(tmp_path):
    source = tmp_path / "benchmark_tools"
    source.mkdir()
    names = ("run_simulation_methods.py", "verify_simulation_histories.py", "run_simulation_generation.py", "benchmark_production.py")
    for name in names:
        (source / name).write_text("pass\n")
    status = {"dataset": "case", "verified_inputs": {"status": "ready"}, "status": "finished_pending_native_validation",
              "provenance": {"method_manifest_sha256": "method", "generation_manifest_sha256": "generation",
                             "slurm_job_id": "11", "slurm_array_task_id": "0",
                             "sources": [file_record(source / n, source) for n in names]}}
    task = {"JobID": "10_0", "JobIDRaw": "11", "State": "COMPLETED"}
    verify_execution_status(status, {"label": "case"}, task, {"status": "ready"}, "method", "generation", tmp_path)
    status["provenance"]["slurm_job_id"] = "12"
    with pytest.raises(ValueError, match="different scheduler"):
        verify_execution_status(status, {"label": "case"}, task, {"status": "ready"}, "method", "generation", tmp_path)


def test_assembler_counts_cross_family_false_positive_and_preserves_failure(tmp_path, monkeypatch):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "s1.fasta").write_text(">a\nAAA\n")
    (inputs / "s2.fasta").write_text(">b\nBBB\n>c\nCCC\n")
    truth = {"extant_genes": 3, "species": ["s1", "s2"], "ortholog_pairs": [["a", "b"]], "ortholog_pair_count": 1}
    truth_path = tmp_path / "truth.json"
    truth_path.write_text(json.dumps(truth))
    verified = {"status": "ready", "truth": file_record(truth_path, tmp_path),
                "inputs": [dict(file_record(p, inputs), absolute_path=str(p)) for p in inputs.iterdir()]}
    methods = {m: {"output": str(tmp_path / "results" / m)} for m in METHODS}
    evidence = tmp_path / "results/execution"
    evidence.mkdir(parents=True)
    artifact = tmp_path / "pairs.txt"
    artifact.write_text("a c\n")
    outputs = [dict(file_record(artifact, tmp_path), absolute_path=str(artifact))]
    (evidence / "status.json").write_text(json.dumps({"methods": {m: {"outputs": outputs} for m in METHODS}}))
    monkeypatch.setattr("benchmark_tools.assemble_simulation_results.verify_inputs", lambda *a: verified)
    monkeypatch.setattr("benchmark_tools.assemble_simulation_results.verify_execution_status", lambda *a: None)
    monkeypatch.setattr("benchmark_tools.assemble_simulation_results.admit_method",
                        lambda m, *a: {"status": "failed", "reason": "native failure"} if m == METHODS[1] else {"status": "admitted"})
    monkeypatch.setattr("benchmark_tools.assemble_simulation_results.load_predictions", lambda *a: ([("a", "c")], [artifact]))
    dataset = {"methods": methods, "truth": str(truth_path), "condition": "baseline", "seed": 20261001}
    rows = dataset_records(dataset, {}, {}, "method", {}, "generation", tmp_path, tmp_path)
    assert len(rows) == 4
    assert rows[0]["score"]["fp"] == 1 and rows[0]["score"]["fn"] == 1 and rows[0]["score"]["tp"] == 0
    assert rows[1]["status"] == "failed" and "score" not in rows[1]
    assert rows[3]["independent_timing"] is False

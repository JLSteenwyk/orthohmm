import json

import pytest

import benchmark_tools.admit_qfo_corrected_factorial_assessment as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("key", ["index", "cell", "stage", "source", "pairs_manifest", "environment_manifest",
    "conversion_scheduler", "converter_commit", "command", "cwd", "verified_records", "environment_overrides"])
def test_reconstructed_provenance_is_required(key):
    expected = {key: "frozen"}
    module.compare_execution(expected, expected)
    with pytest.raises(ValueError, match=key):
        module.compare_execution({key: "changed"}, expected)


def test_reject_nonterminal_before_reading_outputs(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n123|RUNNING|0:0|00:01:00|bizon|8\n")
    with pytest.raises(ValueError, match="COMPLETED"):
        module.admit(tmp_path, 0, "123", "122", "0" * 64, tmp_path / "admission.json")
    assert not (tmp_path / "admission.json").exists()


@pytest.mark.parametrize("index", [0, 1])
def test_admission_orchestration_with_mocked_native_endpoint_checks(tmp_path, monkeypatch, index):
    def write(path, text):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return record(path)
    kind, commit, semantics = module.CONVERTERS[index % 2]
    cell = module.CELLS[index]
    converter = tmp_path / f"benchmarks/work/publication_qfo_corrected_{kind}_pairs_v1"
    converter_source = write(converter / f"benchmark_tools/prepare_qfo_corrected_{kind}_pairs.py", "# converter\n")
    executor = tmp_path / "benchmarks/work/publication_qfo_corrected_factorial_assessment_v1"
    source = write(executor / "benchmark_tools/run_qfo_corrected_factorial_assessment.py", "# runner\n")
    helper = write(executor / "benchmark_tools/run_qfo_recovered_assessment.py", "# command helper\n")
    mapping = write(tmp_path / "mapping.json.gz", "mapping fixture")
    pair = write(tmp_path / "pairs.tsv", "A\tB\n")
    filtered = write(tmp_path / "filtered.tsv", "A\tB\n")
    recheck = write(tmp_path / "native_recheck.json", "{}")
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "8", "JobIDRaw": "123"}
    conversion_scheduler = {**scheduler, "AllocCPUS": "2", "JobIDRaw": "122"}
    stage = {"status": f"corrected_factorial_{kind}_pairs_prepared_unscored", "index": index, "cell": cell,
        "participant": "ohmm_qfo_corrected_factorial_" + cell, "semantics": semantics,
        "accuracy_evaluated": False, "publication_ready": False, "job_id": "122", "total_pairs": 1,
        "retained_pairs": 1, "removed_mapping_pairs": 0, "expected_pairs": 1,
        "pairs": pair, "filtered_pairs": filtered, "mapping": mapping,
        "native_admission_recheck": recheck, "checked_records": [converter_source, recheck]}
    pairs_path = tmp_path / "benchmarks/results/qfo_corrected_factorial_pairs_v1" / cell / "results.json"
    pairs_record = write(pairs_path, json.dumps(stage))
    environment = {"reference_files": [mapping], "pipeline": "/pipeline",
        "execution_config": {"path": "/config"}, "environment_overrides": {"OMP_NUM_THREADS": "1"}}
    env_path = tmp_path / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    env_record = write(env_path, json.dumps(environment))
    monkeypatch.setattr(module, "ENV_SHA", env_record["sha256"])
    monkeypatch.setattr(module, "environment_records", lambda m: [])
    monkeypatch.setattr(module, "verify_checkout", lambda *a: None)
    monkeypatch.setattr(module, "accounting", lambda job: ("assessment", scheduler) if job == "123" else ("conversion", conversion_scheduler))
    results = tmp_path / "qfo_benchmark/scoring" / f"corrected_factorial_{index}"
    trace = write(results / "stats/trace_1.txt", "native trace fixture")
    metric = write(results / "metric.json", "{}")
    directory = tmp_path / "benchmarks/results/qfo_corrected_factorial_assessment_v1" / cell
    log = write(directory / "scoring.log", "completed")
    work = tmp_path / "qfo_benchmark/w" / f"qcf{index}"
    # Isolate orchestration from the production-only Darwin path limit.
    monkeypatch.setattr(module, "command_for", lambda *a: ["nextflow", "--participant_id", stage["participant"]])
    preflight = {"status": "running", "job_id": "123", "accuracy_admitted": False, "index": index,
        "cell": cell, "stage": stage, "source": source, "pairs_manifest": pairs_record,
        "environment_manifest": env_record, "conversion_scheduler": conversion_scheduler,
        "conversion_accounting": "conversion", "converter_commit": commit,
        "command": module.command_for(None), "cwd": str(directory), "work": str(work), "results": str(results),
        "verified_records": [pairs_record, env_record, *stage["checked_records"], pair, filtered, helper],
        "environment_overrides": environment["environment_overrides"]}
    report = {**preflight, "status": "process_succeeded_pending_independent_admission", "exit_code": 0,
              "log": log, "outputs": [metric, trace]}
    write(directory / "preflight.json", json.dumps(preflight))
    write(directory / "results.json", json.dumps(report))
    calls = []
    monkeypatch.setattr(module, "validate_trace", lambda text: calls.append("trace") or ["task"])
    monkeypatch.setattr(module, "validate_directory", lambda path, participant, ref:
        calls.append(participant) or ({"participant": participant}, [results / "metric.json"]))
    output = tmp_path / "admitted.json"
    admitted = module.admit(tmp_path, index, "123", "122", pairs_record["sha256"], output)
    assert admitted["status"] == "corrected_factorial_assessment_admitted"
    assert admitted["accuracy_admitted"] is True and admitted["publication_ready"] is False
    assert calls == ["trace", stage["participant"]]
    with pytest.raises(FileExistsError):
        module.admit(tmp_path, index, "123", "122", pairs_record["sha256"], output)
    (results / "metric.json").write_text("changed")
    with pytest.raises(ValueError, match="output"):
        module.admit(tmp_path, index, "123", "122", pairs_record["sha256"], tmp_path / "changed.json")
    assert not (tmp_path / "changed.json").exists()

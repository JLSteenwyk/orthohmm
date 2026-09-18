import json

import pytest

from benchmark_tools import admit_qfo_corrected_comparator_assessment as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("method", module.METHODS)
def test_method_executor_identity(tmp_path, method):
    executor, commit = module.executor_identity(tmp_path, method)
    if method == "orthomcl":
        assert commit == module.ORTHOMCL_EXECUTOR
        assert executor.name == "publication_qfo_corrected_orthomcl_assessment_v1"
        return
    if method == "fastoma":
        assert commit == module.FASTOMA_EXECUTOR
        assert executor.name == "publication_qfo_corrected_fastoma_assessment_v1"
        return
    is_of = method.startswith("orthofinder_")
    assert commit == (module.OF_EXECUTOR if is_of else module.EXECUTOR)
    assert executor.name == ("publication_qfo_corrected_of_assessment_v1" if is_of
                             else "publication_qfo_corrected_assessment_v1")


@pytest.mark.parametrize("method", ["orthofinder_full", "orthomcl"])
def test_pending_score_refused_before_file_access(tmp_path, monkeypatch, method):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n123|PENDING|0:0|0:00|bizon|8\n")
    with pytest.raises(ValueError, match="COMPLETED"):
        module.admit(tmp_path, method, 123, 122, "0" * 64, tmp_path / "admission.json")
    assert not (tmp_path / "admission.json").exists()


@pytest.mark.parametrize("method", module.METHODS)
def test_orchestration_with_mocked_native_metric_validation(tmp_path, monkeypatch, method):
    def write(path, text):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return record(path)

    is_of = method in module.OF_SEMANTICS
    is_bound = method in module.BOUND_SEMANTICS
    executor, commit = module.executor_identity(tmp_path, method)
    source = write(executor / "benchmark_tools/run_qfo_corrected_comparator_assessment.py", "# runner\n")
    helper = write(executor / "benchmark_tools/run_qfo_recovered_assessment.py", "# command\n")
    pair_helper = write(executor / "benchmark_tools/prepare_qfo_corrected_comparator_pairs.py", "# pairs\n")
    converter = write(tmp_path / "converter.py", "# frozen converter\n")
    pairs = write(tmp_path / "pairs.tsv", "A\tB\n")
    filtered = write(tmp_path / "filtered.tsv", "A\tB\n")
    mapping = write(tmp_path / "mapping.json.gz", "mapping fixture")
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "8", "JobIDRaw": "123"}
    conversion_scheduler = {**scheduler, "AllocCPUS": "2", "JobIDRaw": "122"}
    stage = {"status": "corrected_orthofinder_pairs_prepared_unscored" if is_of
             else "corrected_comparator_pairs_prepared_unscored", "method": method,
             "participant": "qfo_corrected_" + method, "accuracy_evaluated": False,
             "publication_ready": False, "semantics": module.OF_SEMANTICS.get(method, "native"),
             "source": converter, "checked_records": [converter], "pairs": pairs,
             "filtered_pairs": filtered, "mapping": mapping, "job_id": "122",
             "total_pairs": 1, "retained_pairs": 1, "removed_mapping_pairs": 0}
    if method == "fastoma":
        stage.update(status="corrected_fastoma_pairs_prepared_unscored",
                     semantics=module.BOUND_SEMANTICS[method], native_pair_rows=2,
                     native_duplicate_relations=1)
    extra = []
    if method == "orthomcl":
        conversion_scheduler["ReqMem"] = "64G"
        group_audit = write(tmp_path / "groups.json", "audit fixture")
        extra = [group_audit]
        stage.update(status="corrected_orthomcl_pairs_prepared_unscored",
                     semantics=module.BOUND_SEMANTICS[method], native_duplicate_relations=0,
                     content={"total_pairs": 1, "final_groups": 1, "grouped_proteins": 2,
                              "ungrouped_input_proteins": 984135},
                     query_coverage={"failures": ["retained"]}, group_audit=group_audit)
    monkeypatch.setattr(module, "extra_stage_records", lambda *args: extra)
    pairs_path = tmp_path / "benchmarks/results/qfo_corrected_comparator_pairs_v1" / method / "results.json"
    pairs_record = write(pairs_path, json.dumps(stage))
    environment = {"reference_files": [mapping], "pipeline": "/pipeline", "environment_overrides": {"OMP_NUM_THREADS": "1"}}
    env_path = tmp_path / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    env_record = write(env_path, json.dumps(environment))
    monkeypatch.setattr(module, "ENV_SHA", env_record["sha256"])
    monkeypatch.setattr(module, "environment_records", lambda m: [])
    monkeypatch.setattr(module, "converter_source", lambda root, method: converter if is_bound else None)
    def fake_accounting(job, include_memory=False):
        if job != "123":
            assert include_memory == (method == "orthomcl")
        return ("score", scheduler) if job == "123" else ("conversion", conversion_scheduler)
    monkeypatch.setattr(module, "accounting", fake_accounting)
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: commit + "\n")
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    directory = tmp_path / "benchmarks/results/qfo_corrected_assessment_v1" / method
    results = tmp_path / "qfo_benchmark/scoring" / ("corrected_" + method)
    work = tmp_path / "qfo_benchmark/w" / module.WORK_NAMES[method]
    metric = write(results / "metric.json", "{}")
    trace = write(results / "stats/trace_1.txt", "native trace fixture")
    log = write(directory / "scoring.log", "completed")
    monkeypatch.setattr(module, "command_for", lambda *a: ["nextflow", stage["participant"]])
    checked = [pairs_record, env_record, converter, converter, pairs, filtered, helper, pair_helper]
    if is_bound:
        checked.append(converter)
    checked.extend(extra)
    preflight = {"status": "running", "job_id": "123", "accuracy_admitted": False,
        "method": method, "stage": stage, "source": source, "pairs_manifest": pairs_record,
        "environment_manifest": env_record, "conversion_scheduler": conversion_scheduler,
        "conversion_accounting": "conversion", "command": module.command_for(None),
        "cwd": str(directory), "work": str(work), "results": str(results),
        "verified_records": checked, "environment_overrides": environment["environment_overrides"]}
    report = {**preflight, "status": "process_succeeded_pending_independent_admission",
              "exit_code": 0, "log": log, "outputs": [metric, trace]}
    write(directory / "preflight.json", json.dumps(preflight))
    write(directory / "results.json", json.dumps(report))
    calls = []
    monkeypatch.setattr(module, "validate_trace", lambda text: calls.append("trace") or ["task"])
    monkeypatch.setattr(module, "validate_directory", lambda path, participant, ref:
        calls.append(participant) or ({"participant": participant}, [results / "metric.json"]))
    output = tmp_path / "admitted.json"
    admitted = module.admit(tmp_path, method, "123", "122", pairs_record["sha256"], output)
    assert admitted["accuracy_admitted"] is True and admitted["publication_ready"] is False
    assert calls == ["trace", stage["participant"]]
    assert admitted["execution_report"] == record(directory / "results.json")
    if method == "orthomcl":
        assert admitted["query_coverage"] == stage["query_coverage"]
        assert admitted["group_coverage"] == stage["content"]
        assert admitted["pair_semantics"] == module.BOUND_SEMANTICS[method]
        assert group_audit in admitted["checked_records"]
    with pytest.raises(FileExistsError):
        module.admit(tmp_path, method, "123", "122", pairs_record["sha256"], output)
    (results / "metric.json").write_text("changed")
    with pytest.raises(ValueError, match="output"):
        module.admit(tmp_path, method, "123", "122", pairs_record["sha256"], tmp_path / "changed.json")
    assert not (tmp_path / "changed.json").exists()
    if is_bound:
        environment["reference_files"] = []
        changed_env = write(env_path, json.dumps(environment))
        monkeypatch.setattr(module, "ENV_SHA", changed_env["sha256"])
        with pytest.raises(ValueError, match="reference mappings differ"):
            module.admit(tmp_path, method, "123", "122", pairs_record["sha256"],
                         tmp_path / "wrong_mapping.json")
        assert not (tmp_path / "wrong_mapping.json").exists()

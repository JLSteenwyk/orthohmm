import pytest

from benchmark_tools.run_qfo_corrected_comparator_assessment import validate_stage


def fixture():
    stage = {"status": "corrected_comparator_pairs_prepared_unscored", "accuracy_evaluated": False,
             "method": "proteinortho", "participant": "qfo_corrected_proteinortho", "job_id": "1",
             "total_pairs": 2, "retained_pairs": 2, "removed_mapping_pairs": 0,
             "pairs": {"bytes": 8, "sha256": "a"}, "filtered_pairs": {"bytes": 8, "sha256": "a"}}
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "2", "JobIDRaw": "1"}
    return stage, scheduler


def test_valid():
    stage, scheduler = fixture()
    validate_stage(stage, "proteinortho", scheduler)


@pytest.mark.parametrize("key,value", [("status", "failed"), ("accuracy_evaluated", True),
    ("method", "sonic"), ("participant", "historical"), ("job_id", "2"),
    ("total_pairs", 0), ("total_pairs", True), ("retained_pairs", 1), ("removed_mapping_pairs", 1),
    ("filtered_pairs", {"bytes": 8, "sha256": "changed"})])
def test_reject_changed_stage(key, value):
    stage, scheduler = fixture()
    stage[key] = value
    with pytest.raises(ValueError):
        validate_stage(stage, "proteinortho", scheduler)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
                                      ("NodeList", "other"), ("AllocCPUS", "1")])
def test_reject_scheduler(key, value):
    stage, scheduler = fixture()
    scheduler[key] = value
    with pytest.raises(ValueError, match="scheduler"):
        validate_stage(stage, "proteinortho", scheduler)


@pytest.mark.parametrize("exit_code", [0, 1])
def test_run_preserves_native_status(monkeypatch, tmp_path, exit_code):
    import json
    from types import SimpleNamespace
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    source = tmp_path / "source"
    source.write_text("fixture")
    output = tmp_path / "execution"
    report = {"cwd": str(output), "results": str(tmp_path / "scores"), "command": ["fixture"],
              "source": module.record(source), "verified_records": [], "environment_overrides": {},
              "accuracy_admitted": False}
    monkeypatch.setattr(module, "prepare", lambda *args: report)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=exit_code))
    if exit_code:
        with pytest.raises(RuntimeError, match="scoring failed"):
            module.run(tmp_path, "proteinortho", "fixture", 1)
    else:
        assert module.run(tmp_path, "proteinortho", "fixture", 1)["status"] == "process_succeeded_pending_independent_admission"
    final = json.loads((output / "results.json").read_text())
    assert final["exit_code"] == exit_code
    assert final["accuracy_admitted"] is False
    assert (output / "preflight.json").exists()
    assert final["status"] == ("failed" if exit_code else "process_succeeded_pending_independent_admission")


def orthofinder_fixture(method):
    from benchmark_tools.run_qfo_corrected_comparator_assessment import OF_SEMANTICS
    stage, scheduler = fixture()
    stage.update(method=method, participant="qfo_corrected_" + method,
        status="corrected_orthofinder_pairs_prepared_unscored", publication_ready=False,
        semantics=OF_SEMANTICS[method])
    return stage, scheduler


@pytest.mark.parametrize("method", ["orthofinder_full", "orthofinder_sequence_only"])
def test_orthofinder_semantics(method):
    stage, scheduler = orthofinder_fixture(method)
    validate_stage(stage, method, scheduler)


@pytest.mark.parametrize("key,value", [("status", "corrected_comparator_pairs_prepared_unscored"),
    ("semantics", "root HOG cliques"), ("publication_ready", True),
    ("method", "orthofinder_sequence_only"), ("participant", "historical_of"),
    ("retained_pairs", True), ("removed_mapping_pairs", False)])
def test_wrong_orthofinder_stage_refused(key, value):
    stage, scheduler = orthofinder_fixture("orthofinder_full")
    stage[key] = value
    with pytest.raises(ValueError):
        validate_stage(stage, "orthofinder_full", scheduler)


@pytest.mark.parametrize("method", ["orthofinder_full", "fastoma", "orthomcl"])
def test_pending_conversion_rejected_before_reading_partial_output(tmp_path, monkeypatch, method):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n123|PENDING|0:0|0:00|bizon|2\n")
    monkeypatch.setattr(module, "read_frozen", lambda *a: pytest.fail("Partial conversion read"))
    with pytest.raises(ValueError, match="COMPLETED"):
        module.prepare(tmp_path, method, "sha", 123)


def test_all_workspaces_are_distinct():
    from benchmark_tools.run_qfo_corrected_comparator_assessment import METHODS, WORK_NAMES
    assert len(set(WORK_NAMES.values())) == len(METHODS) == 6
    assert WORK_NAMES["proteinortho"] == "qc_p"
    assert WORK_NAMES["sonic"] == "qc_s"
    assert WORK_NAMES["orthomcl"] == "qc_mc"


def orthomcl_fixture():
    from benchmark_tools.run_qfo_corrected_comparator_assessment import ORTHOMCL_SEMANTICS
    stage, scheduler = fixture()
    stage.update(status="corrected_orthomcl_pairs_prepared_unscored", method="orthomcl",
                 participant="qfo_corrected_orthomcl", publication_ready=False,
                 semantics=ORTHOMCL_SEMANTICS, native_duplicate_relations=0,
                 content={"total_pairs": 2, "final_groups": 1, "grouped_proteins": 3,
                          "ungrouped_input_proteins": 984134})
    scheduler["ReqMem"] = "64G"
    return stage, scheduler


def test_orthomcl_stage():
    stage, scheduler = orthomcl_fixture()
    validate_stage(stage, "orthomcl", scheduler)


@pytest.mark.parametrize("key,value", [("status", "corrected_comparator_pairs_prepared_unscored"),
    ("semantics", "pre-clustering graph edges"), ("publication_ready", True),
    ("native_duplicate_relations", 1), ("native_duplicate_relations", False), ("content", {}),
    ("participant", "historical_orthomcl"), ("removed_mapping_pairs", 1)])
def test_bad_orthomcl_stage(key, value):
    stage, scheduler = orthomcl_fixture()
    stage[key] = value
    with pytest.raises(ValueError):
        validate_stage(stage, "orthomcl", scheduler)


@pytest.mark.parametrize("key,value", [("total_pairs", 3), ("total_pairs", True),
    ("final_groups", 0), ("grouped_proteins", 984138), ("ungrouped_input_proteins", -1)])
def test_bad_orthomcl_counts(key, value):
    stage, scheduler = orthomcl_fixture()
    stage["content"][key] = value
    with pytest.raises(ValueError):
        validate_stage(stage, "orthomcl", scheduler)


def test_orthomcl_conversion_memory_bound():
    stage, scheduler = orthomcl_fixture()
    scheduler["ReqMem"] = "32G"
    with pytest.raises(ValueError):
        validate_stage(stage, "orthomcl", scheduler)


def test_orthomcl_converter_revision_pinned(tmp_path, monkeypatch):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "wrong\n")
    with pytest.raises(ValueError, match="executor changed"):
        module.converter_source(tmp_path, "orthomcl")


def test_orthomcl_converter_source(tmp_path, monkeypatch):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    executor = tmp_path / "benchmarks/work/publication_qfo_corrected_orthomcl_pairs_v1"
    source = executor / "benchmark_tools/prepare_qfo_corrected_orthomcl_pairs.py"
    source.parent.mkdir(parents=True)
    source.write_text("fixture")
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: module.ORTHOMCL_CONVERTER + "\n")
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    assert module.converter_source(tmp_path, "orthomcl") == module.record(source)


def orthomcl_evidence(tmp_path):
    import json
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    stage, scheduler = orthomcl_fixture()
    directory = tmp_path / "benchmarks/results/qfo_corrected_comparator_pairs_v1/orthomcl"
    directory.mkdir(parents=True)
    content = {"input_proteins": 984137, "input_species": 78, "cross_species_clique_pairs": 2,
               **{k: v for k, v in stage["content"].items() if k != "total_pairs"}}
    native = {"status": "corrected_orthomcl_native_outputs_admitted", "accuracy_admitted": False,
              "publication_ready": False, "pair_semantics": module.ORTHOMCL_SEMANTICS,
              "content": content, "query_coverage": {"failures": ["retained"]},
              "scheduler": {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "180", "ReqMem": "900G"}}
    path = tmp_path / "benchmarks/work/qfo_corrected_orthomcl_admission_20260918/report.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(native))
    groups = {"status": "native_final_groups_match_mcl_partition", "content": content,
              "accuracy_admitted": False, "publication_ready": False, "checked_records": []}
    (directory / "groups.json").write_text(json.dumps(groups))
    for name in ("pairs.tsv", "pairs.qfo.tsv"):
        (directory / name).write_text("A\tB\nA\tC\n")
    for key, path in (("admission", path), ("group_audit", directory / "groups.json"),
                      ("pairs", directory / "pairs.tsv"), ("filtered_pairs", directory / "pairs.qfo.tsv")):
        stage[key] = module.record(path)
    stage["checked_records"] = [stage["admission"]]
    stage["query_coverage"] = native["query_coverage"]
    return stage


def test_orthomcl_additional_evidence(tmp_path):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    stage = orthomcl_evidence(tmp_path)
    assert module.extra_stage_records(tmp_path, "orthomcl", stage) == [stage["group_audit"], stage["admission"]]
    assert module.extra_stage_records(tmp_path, "fastoma", {}) == []


@pytest.mark.parametrize("change", ["path", "unbound", "diagnostics", "count", "audit", "source_hash"])
def test_orthomcl_extra_evidence_rejects_changes(tmp_path, change):
    import json
    from pathlib import Path
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    stage = orthomcl_evidence(tmp_path)
    if change == "path":
        stage["pairs"]["path"] = str(tmp_path / "unbound.tsv")
    elif change == "unbound":
        stage["checked_records"] = []
    elif change == "diagnostics":
        stage["query_coverage"] = {}
    elif change == "count":
        stage["content"]["final_groups"] += 1
    elif change == "source_hash":
        Path(stage["admission"]["path"]).write_text("{}")
    else:
        path = Path(stage["group_audit"]["path"])
        data = json.loads(path.read_text())
        data["content"]["cross_species_clique_pairs"] += 1
        path.write_text(json.dumps(data))
        stage["group_audit"] = module.record(path)
    with pytest.raises(ValueError):
        module.extra_stage_records(tmp_path, "orthomcl", stage)


def test_unknown_method_rejected(tmp_path):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    with pytest.raises(ValueError, match="Unknown"):
        module.prepare(tmp_path, "unknown", "sha", 123)


def test_converter_revision_pinned(tmp_path, monkeypatch):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "wrong\n")
    with pytest.raises(ValueError, match="executor changed"):
        module.converter_source(tmp_path, "orthofinder_full")


def fastoma_fixture():
    from benchmark_tools.run_qfo_corrected_comparator_assessment import FASTOMA_SEMANTICS
    stage, scheduler = fixture()
    stage.update(status="corrected_fastoma_pairs_prepared_unscored", method="fastoma",
                 participant="qfo_corrected_fastoma", publication_ready=False,
                 semantics=FASTOMA_SEMANTICS, native_pair_rows=3, native_duplicate_relations=1)
    return stage, scheduler


def test_fastoma_stage():
    stage, scheduler = fastoma_fixture()
    validate_stage(stage, "fastoma", scheduler)


@pytest.mark.parametrize("key,value", [("status", "corrected_comparator_pairs_prepared_unscored"),
    ("semantics", "HOG cliques"), ("publication_ready", True), ("native_pair_rows", 2),
    ("native_duplicate_relations", -1), ("native_duplicate_relations", True),
    ("participant", "old_fastoma"), ("removed_mapping_pairs", 1)])
def test_bad_fastoma_stage(key, value):
    stage, scheduler = fastoma_fixture()
    stage[key] = value
    with pytest.raises(ValueError):
        validate_stage(stage, "fastoma", scheduler)


def test_fastoma_converter_revision_pinned(tmp_path, monkeypatch):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "wrong\n")
    with pytest.raises(ValueError, match="executor changed"):
        module.converter_source(tmp_path, "fastoma")


def test_fastoma_converter_source_is_frozen_worktree(tmp_path, monkeypatch):
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    executor = tmp_path / "benchmarks/work/publication_qfo_corrected_fastoma_pairs_v2"
    source = executor / "benchmark_tools/prepare_qfo_corrected_fastoma_pairs.py"
    source.parent.mkdir(parents=True)
    source.write_text("fixture")
    calls = []
    monkeypatch.setattr(module.subprocess, "check_output",
                        lambda *a, **k: module.FASTOMA_CONVERTER + "\n")
    monkeypatch.setattr(module.subprocess, "run", lambda args, **k: calls.append((args, k)))
    assert module.converter_source(tmp_path, "fastoma") == module.record(source)
    assert calls == [(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--",
                      "benchmark_tools"], {"check": True})]

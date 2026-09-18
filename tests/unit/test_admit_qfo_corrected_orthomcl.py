import copy
import json
import os
from pathlib import Path

import pytest

from benchmark_tools import admit_qfo_corrected_orthomcl as module


def execution(tmp_path):
    executor, base = tmp_path / "executor", tmp_path / "base"
    source = executor / "benchmark_tools/run_qfo_corrected_orthomcl.py"
    source.parent.mkdir(parents=True)
    source.write_text("source\n")
    directory = base / "inference_execution"
    (directory / "inputs").mkdir(parents=True)
    for name in ("native.log", "native.time.txt", "inputs/staging.json"):
        (directory / name).write_text("fixture\n")
    scheduler = {"JobIDRaw": "123", "State": "COMPLETED", "ExitCode": "0:0",
                 "NodeList": "bizon", "AllocCPUS": "180", "ReqMem": "900G"}
    runtime = {"path": "pinned", "sha256": "abc", "bytes": 1}
    tool = base / "native_tool"
    report = {"status": "corrected_orthomcl_native_exited_zero_pending_admission", "job_id": "123",
              "node": "bizon", "native_exit_code": 0, "accuracy_admitted": False, "publication_ready": False,
              "source": module.record(source), "cwd": str(directory),
              "environment": {**module.environment(), "ORTHOMCL_PAIR_WORKERS": "64"},
              "started_epoch": 1, "native_started_epoch": 2, "native_finished_epoch": 3, "finished_epoch": 4,
              "command": [str(module.TOOL / "venv_orthomcl/bin/perl"), "-I" + str(tool),
                  str(executor / "benchmark_tools/run_orthomcl_perl_script.pl"), str(tool / "orthomcl.pl"),
                  "--mode", "4", "--bpo_file", str(directory / "inputs/all.bpo"),
                  "--gg_file", str(directory / "inputs/all.gg")]}
    for key, name in (("native_log", "native.log"), ("native_timing", "native.time.txt"), ("staging", "inputs/staging.json")):
        report[key] = module.record(directory / name)
    for key in ("runtime_before", "runtime_after"):
        report[key] = {"status": "dedicated_bpo_python_runtime_verified", "manifest": runtime, "mapped_files": []}
    return report, scheduler, executor, base, runtime


def test_execution_valid(tmp_path):
    module.validate_execution(*execution(tmp_path))


@pytest.mark.parametrize("key,value", [
    ("status", "failed"), ("job_id", "124"), ("node", "spark"), ("native_exit_code", 1),
    ("native_exit_code", False), ("accuracy_admitted", True), ("publication_ready", True),
    ("cwd", "/elsewhere"), ("command", []), ("environment", {}), ("source", {}),
    ("native_log", {}), ("native_timing", {}), ("staging", {}),
    ("started_epoch", 0), ("finished_epoch", 1), ("native_started_epoch", 5),
    ("native_finished_epoch", float("nan")), ("finished_epoch", float("inf")),
    ("started_epoch", True),
])
def test_execution_corruption(tmp_path, key, value):
    args = execution(tmp_path)
    args[0][key] = value
    with pytest.raises(ValueError):
        module.validate_execution(*args)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"), ("NodeList", "spark"),
                                      ("AllocCPUS", "64"), ("ReqMem", "64G")])
def test_scheduler_corruption(tmp_path, key, value):
    args = execution(tmp_path)
    args[1][key] = value
    with pytest.raises(ValueError):
        module.validate_execution(*args)


@pytest.mark.parametrize("stage", ["runtime_before", "runtime_after"])
def test_runtime_manifest_corruption(tmp_path, stage):
    args = execution(tmp_path)
    args[0][stage]["manifest"] = {}
    with pytest.raises(ValueError, match="runtime"):
        module.validate_execution(*args)


def staging_fixture(tmp_path):
    directory = tmp_path / "inputs"
    (directory / "validation").mkdir(parents=True)
    executor = tmp_path / "executor"
    inputs, staged, mtimes = {}, {}, {}
    for key, name in module.NAMES.items():
        source = tmp_path / (name + ".source")
        source.write_text(key + "\n")
        target = directory / name
        target.write_bytes(source.read_bytes())
        inputs[key], staged[key] = module.record(source), module.record(target)
        mtimes[key] = target.stat().st_mtime_ns
    indexes = {"records": 4, "queries": 2}
    (directory / "validation/indexes.stdout").write_text(json.dumps(indexes))
    (directory / "validation/indexes.stderr").write_text("")
    staging = {"status": "native_inputs_staged_and_indexes_verified", "accuracy_admitted": False,
               "publication_ready": False, "sources": inputs, "staged": staged,
               "input_proteins": 984137, "species": 78, "index_validation": indexes,
               "index_validation_command": [str(module.TOOL / "venv_orthomcl/bin/perl"),
                   str(executor / "benchmark_tools/run_orthomcl_perl_script.pl"),
                   str(executor / "benchmark_tools/validate_orthomcl_bpo_indexes.pl"),
                   *[str(directory / module.NAMES[k]) for k in ("bpo", "offsets", "ranges")]]}
    report = {"staged_inputs": staged, "input_mtimes_ns": mtimes}
    return staging, report, inputs, indexes, directory, executor


def test_staging_valid(tmp_path):
    module.validate_staging(*staging_fixture(tmp_path))


@pytest.mark.parametrize("key,value", [("status", "failed"), ("accuracy_admitted", True), ("publication_ready", True),
    ("sources", {}), ("staged", {}), ("input_proteins", 42), ("species", 3), ("index_validation", {}),
    ("index_validation_command", [])])
def test_staging_corruption(tmp_path, key, value):
    args = staging_fixture(tmp_path)
    args[0][key] = value
    with pytest.raises(ValueError):
        module.validate_staging(*args)


@pytest.mark.parametrize("change", ["mtime", "bytes", "hardlink", "symlink", "stderr", "stdout"])
def test_changed_staged_files_rejected(tmp_path, change):
    args = staging_fixture(tmp_path)
    path = args[4] / "all.bpo"
    if change == "mtime":
        os.utime(path, ns=(1, 1))
    elif change == "bytes":
        old = path.stat()
        path.write_text("changed")
        os.utime(path, ns=(old.st_atime_ns, old.st_mtime_ns))
    elif change == "hardlink":
        os.link(path, tmp_path / "link")
    elif change == "symlink":
        path.unlink()
        path.symlink_to(args[2]["bpo"]["path"])
    else:
        (args[4] / ("validation/indexes." + change)).write_text("{}")
    with pytest.raises(ValueError):
        module.validate_staging(*args)


def outputs_fixture(tmp_path, monkeypatch):
    base = tmp_path / "base"
    directory = base / "native_tool/Sep_18"
    (directory / "tmp").mkdir(parents=True)
    for name in ("all_orthomcl.out", "tmp/all_ortho.idx", "tmp/all_ortho.mtx", "tmp/all_ortho.mcl"):
        (directory / name).write_text(name + "\n")
    report = {"native_groups": module.record(directory / "all_orthomcl.out"),
              "native_outputs": [module.record(p) for p in sorted(directory.rglob("*")) if p.is_file()],
              "pair_caches": [{"cache": "bound"}]}
    monkeypatch.setattr(module, "pair_cache_records", lambda path: [{"cache": "bound"}])
    return report, base, directory


def test_output_inventory(tmp_path, monkeypatch):
    report, base, directory = outputs_fixture(tmp_path, monkeypatch)
    assert module.validate_outputs(report, base) == directory


@pytest.mark.parametrize("change", ["extra", "missing", "duplicate", "symlink", "cache", "modified"])
def test_bad_output_inventory(tmp_path, monkeypatch, change):
    report, base, directory = outputs_fixture(tmp_path, monkeypatch)
    if change == "extra":
        (directory / "extra").write_text("extra")
    elif change == "missing":
        (directory / "tmp/all_ortho.mcl").unlink()
    elif change == "duplicate":
        other = base / "native_tool/other"
        other.mkdir()
        (other / "all_orthomcl.out").write_text("duplicate")
    elif change == "symlink":
        (directory / "link").symlink_to(directory / "all_orthomcl.out")
    elif change == "cache":
        report["pair_caches"] = []
    else:
        (directory / "all_orthomcl.out").write_text("changed")
    with pytest.raises(ValueError):
        module.validate_outputs(report, base)


def test_conflicting_records_rejected():
    item = {"path": "x", "sha256": "one", "bytes": 1}
    assert module.unique_records([item, copy.deepcopy(item)]) == [item]
    with pytest.raises(ValueError, match="Conflicting"):
        module.unique_records([item, {**item, "sha256": "two"}])


def test_pending_job_fails_before_reading_outputs(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "verify_runtime", lambda root: {})
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n123|PENDING|0:0|00:00:00||180|900G\n")
    output = tmp_path / "admission"
    with pytest.raises(ValueError):
        module.admit(tmp_path, 123, output)
    assert not output.exists()


def test_no_overwrite_dangling_symlink(tmp_path):
    output = tmp_path / "report"
    output.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        module.admit(tmp_path, 123, output)


@pytest.mark.parametrize("failure", [None, "groups", "indexes"])
def test_admission_orchestration_and_retained_failures(tmp_path, monkeypatch, failure):
    root = tmp_path
    base = root / "benchmarks/results/qfo_corrected_orthomcl_v1"
    directory = base / "inference_execution/inputs"
    directory.mkdir(parents=True)
    native = base / "native_tool/Sep_18"
    (native / "tmp").mkdir(parents=True)
    for path in (native / "all_orthomcl.out", directory / "all.gg", directory / "all.bpo",
                 directory / "all_bpo.idx", directory / "all_bpo.se"):
        path.write_text("fixture")
    inputs = {key: module.record(directory / name) for key, name in module.NAMES.items()}
    stage_path = directory / "staging.json"
    stage_path.write_text("{}")
    executor = root / "benchmarks/work/publication_qfo_corrected_orthomcl_v1"
    prior = root / "benchmarks/work/publication_qfo_corrected_bpo_admission_v1"
    for worktree, name in ((executor, "run_qfo_corrected_orthomcl.py"), (prior, "admit_qfo_corrected_bpo.py")):
        source = worktree / "benchmark_tools" / name
        source.parent.mkdir(parents=True)
        source.write_text("fixture")
    prior_scheduler = {"JobIDRaw": "122"}
    admission = {"source": module.record(prior / "benchmark_tools/admit_qfo_corrected_bpo.py"),
                 "query_coverage": {"failed_query_ids": ["retained-failure"]}, "checked_records": [],
                 "validation": {"outputs": [], "index_validation": {"queries": 3}}}
    admission_path = root / "benchmarks/work/qfo_corrected_bpo_admission_20260918/report.json"
    admission_path.parent.mkdir(parents=True)
    admission_path.write_text(json.dumps(admission))
    results = root / "benchmark_tools/results"
    results.mkdir(parents=True)
    source_path = results / "qfo_corrected_orthomcl_native_sources_20260918.json"
    source_path.write_text(json.dumps({"checked_records": [], "originals": [], "configured_sources": []}))
    for name in ("qfo_corrected_orthomcl_perl_runtime_20260918.json", "qfo_corrected_orthomcl_system_helpers_20260918.json"):
        (results / name).write_text("{}")
    report = {"checked_records": [module.record(admission_path)], "staging": module.record(stage_path),
              "admission_scheduler": prior_scheduler, "query_coverage": admission["query_coverage"],
              "native_log": module.record(stage_path), "native_timing": module.record(stage_path),
              "staged_inputs": inputs, "native_groups": module.record(native / "all_orthomcl.out"),
              "content": {"groups": 2, "grouped_proteins": 5, "ungrouped_proteins": 984132}}
    (base / "inference_execution/status.json").write_text(json.dumps(report))
    content = {"input_proteins": 984137, "input_species": 78, "final_groups": 2,
               "grouped_proteins": 5, "ungrouped_input_proteins": 984132}
    monkeypatch.setattr(module, "verify_runtime", lambda root: {"manifest": {}})
    monkeypatch.setattr(module, "completed", lambda job, *args: (prior_scheduler if job == 122 else {"JobIDRaw": "123"}, "accounting"))
    monkeypatch.setattr(module, "frozen", lambda *args: None)
    monkeypatch.setattr(module, "validate_execution", lambda *args: None)
    monkeypatch.setattr(module, "validate_staging", lambda *args: None)
    monkeypatch.setattr(module, "validate_outputs", lambda *args: native)
    monkeypatch.setattr(module, "admitted_inputs", lambda *args: inputs)
    monkeypatch.setattr(module, "read_frozen", lambda path, sha: json.loads(path.read_text()))
    monkeypatch.setattr(module, "verify", lambda *args: None)
    def fake_audit(*args):
        observed = {**content, "final_groups": 3} if failure == "groups" else content
        result = {"content": observed}
        args[-1].write_text(json.dumps(result))
        return result
    monkeypatch.setattr(module, "audit", fake_audit)
    def fake_indexes(argv, output, env, stage):
        (output / "indexes.stdout").write_text(json.dumps({"queries": 0 if failure == "indexes" else 3}))
        (output / "indexes.stderr").write_text("")
    monkeypatch.setattr(module, "native_step", fake_indexes)
    output = root / "admission"
    if failure:
        with pytest.raises(ValueError):
            module.admit(root, 123, output)
        saved = json.loads((output / "report.json").read_text())
        assert saved["status"] == "failed"
    else:
        saved = module.admit(root, 123, output)
        assert saved["status"] == "corrected_orthomcl_native_outputs_admitted"
        assert saved["pair_semantics"] == "cross_species_final_group_cliques"
        assert saved == json.loads((output / "report.json").read_text())
        assert len(saved["outputs"]) == 3
    assert saved["query_coverage"] == admission["query_coverage"]
    assert saved["accuracy_admitted"] is saved["publication_ready"] is False

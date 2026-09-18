import copy
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from benchmark_tools import run_qfo_corrected_orthomcl as module


def admission(root):
    base = root / "benchmarks/results/qfo_corrected_orthomcl_v1"
    paths = [base / "bpo_preparation/checkpoint" / name for name in ("all.bpo", "indexes/all_bpo.idx", "indexes/all_bpo.se")]
    records = [{"path": str(p), "bytes": 1, "sha256": "x"} for p in [*paths, base / "work/all.gg"]]
    return {"status": "corrected_orthomcl_bpo_checkpoint_admitted", "accuracy_admitted": False,
            "publication_ready": False, "validation": {"content": {"input_proteins": 984137}},
            "scheduler": {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "2", "ReqMem": "64G"},
            "native_inputs": records[:3], "checked_records": records}


def test_admitted_inputs(tmp_path):
    report = admission(tmp_path)
    result = module.admitted_inputs(report, tmp_path)
    assert set(result) == {"bpo", "offsets", "ranges", "species"}
    report["checked_records"].append(copy.deepcopy(report["checked_records"][-1]))
    assert module.admitted_inputs(report, tmp_path) == result


@pytest.mark.parametrize("problem", ["status", "accuracy", "publication", "scope", "state", "allocation", "inventory", "missing_gg", "conflict", "empty", "changed_hash"])
def test_bad_admission(tmp_path, problem):
    report = admission(tmp_path)
    if problem == "status":
        report["status"] = "running"
    elif problem in {"accuracy", "publication"}:
        report["accuracy_admitted" if problem == "accuracy" else "publication_ready"] = True
    elif problem == "scope":
        report["validation"]["content"]["input_proteins"] -= 1
    elif problem in {"state", "allocation"}:
        report["scheduler"]["State" if problem == "state" else "AllocCPUS"] = "wrong"
    elif problem == "inventory":
        report["native_inputs"].append(report["native_inputs"][0])
    elif problem == "missing_gg":
        report["checked_records"].pop()
    elif problem == "conflict":
        report["checked_records"].append({**report["checked_records"][0], "sha256": "bad"})
    elif problem == "empty":
        report["checked_records"][-1]["bytes"] = 0
    else:
        report["native_inputs"] = copy.deepcopy(report["native_inputs"])
        report["native_inputs"][0]["sha256"] = "bad"
    with pytest.raises(ValueError):
        module.admitted_inputs(report, tmp_path)


@pytest.mark.parametrize("text", ["", "OG0: A(a) A(a)\n", "OG0: X(a)\n", "OG0: A(a)\nOG0: B(b)\n", "OG0: A(a)\nOG1: A(a)\n"])
def test_invalid_groups(tmp_path, text):
    path = tmp_path / "groups"
    path.write_text(text)
    with pytest.raises(ValueError):
        module.group_counts(path, {"A": "a", "B": "b"})


def test_counts_preserve_ungrouped_input(tmp_path):
    path = tmp_path / "groups"
    path.write_text("OG0: A(a)\n")
    assert module.group_counts(path, {"A": "a", "B": "b"}) == {"groups": 1, "grouped_proteins": 1, "ungrouped_proteins": 1}


@pytest.mark.parametrize("problem", [None, "missing", "extra", "empty", "symlink"])
def test_pair_cache_inventory(tmp_path, problem):
    cache = tmp_path / "pair_0.storable"
    if problem != "missing":
        cache.write_text("sample" if problem != "empty" else "")
    if problem == "extra":
        (tmp_path / "pair_0.storable.partial").write_text("partial")
    elif problem == "symlink":
        cache.rename(tmp_path / "real")
        cache.symlink_to(tmp_path / "real")
    if problem:
        with pytest.raises(ValueError):
            module.pair_cache_records(tmp_path, 2)
    else:
        assert module.pair_cache_records(tmp_path, 2) == [module.record(cache)]


def test_command_uses_mode4_and_guarded_perl(tmp_path):
    argv = module.command(tmp_path / "tool", tmp_path / "inputs")
    assert argv[2].endswith("run_orthomcl_perl_script.pl")
    assert argv[4:] == ["--mode", "4", "--bpo_file", str(tmp_path / "inputs/all.bpo"), "--gg_file", str(tmp_path / "inputs/all.gg")]


@pytest.mark.parametrize("problem", [None, "native_failure", "missing_groups", "modified_input", "stage_failure"])
def test_execution_state_and_failure_preservation(tmp_path, monkeypatch, problem):
    for key, value in {"SLURM_JOB_ID": "123", "SLURM_CPUS_PER_TASK": "180", "SLURM_MEM_PER_NODE": "921600"}.items():
        monkeypatch.setenv(key, value)
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    tool = tmp_path / "native_tool"
    tool.mkdir()
    monkeypatch.setattr(module, "preflight", lambda *args: ({"query_coverage": {}, "validation": {"index_validation": {}}}, {}, tool, [], [], {}, {}, ""))
    monkeypatch.setattr(module, "verify_runtime", lambda root: {})
    monkeypatch.setattr(module, "pair_cache_records", lambda path: [])
    def stage(inputs, directory, *args):
        directory.mkdir()
        if problem == "stage_failure":
            raise ValueError("staging failed")
        (directory / "staging.json").write_text("{}")
        (directory / "all.gg").write_text("a: A\nb: B\n")
        return {"staged": {"species": module.record(directory / "all.gg")}}
    monkeypatch.setattr(module, "stage", stage)
    def execute(argv, **kwargs):
        assert kwargs["env"]["ORTHOMCL_PAIR_WORKERS"] == "64"
        Path(argv[3]).write_text("timing")
        if problem == "modified_input":
            (tmp_path / "inference_execution/inputs/all.gg").write_text("changed")
        if problem != "missing_groups":
            output = tool / "native"
            output.mkdir()
            (output / "all_orthomcl.out").write_text("OG0: A(a) B(b)\n")
        return subprocess.CompletedProcess(argv, 1 if problem == "native_failure" else 0)
    monkeypatch.setattr(module.subprocess, "run", execute)
    if problem:
        with pytest.raises(ValueError):
            module.run(tmp_path, tmp_path / "admission", "sha", 99)
        report = json.loads((tmp_path / "inference_execution/status.json").read_text())
        assert report["status"] == "failed"
    else:
        report = module.run(tmp_path, tmp_path / "admission", "sha", 99)
        assert report["status"] == "corrected_orthomcl_native_exited_zero_pending_admission"
        assert report["content"] == {"groups": 1, "grouped_proteins": 2, "ungrouped_proteins": 0}
    assert report["accuracy_admitted"] is report["publication_ready"] is False


def test_unscheduled_execution_rejected(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        module.run(tmp_path, tmp_path / "admission", "sha", 99)


@pytest.mark.parametrize("problem", [None, "commit", "admission_source", "allocation", "workers", "old_tool_output", "old_execution", "space", "changed_input"])
def test_preflight_provenance_and_freshness(tmp_path, monkeypatch, problem):
    data = admission(tmp_path)
    for item in data["checked_records"]:
        path = Path(item["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("input")
        item.update(module.record(path))
    data["validation"]["outputs"] = []
    executor = tmp_path / "benchmarks/work/publication_qfo_corrected_bpo_admission_v1"
    (executor / "benchmark_tools").mkdir(parents=True)
    source = executor / "benchmark_tools/admit_qfo_corrected_bpo.py"
    source.write_text("source")
    data["source"] = module.record(source)
    if problem == "admission_source":
        data["source"]["sha256"] = "wrong"
    admission_path = tmp_path / "admission.json"
    admission_path.write_text("{}")
    tool = tmp_path / "benchmarks/results/qfo_corrected_orthomcl_v1/native_tool"
    tool.mkdir()
    for name in ("orthomcl.pl", "orthomcl_module.pm"):
        (tool / name).write_text("native source")
    source_path = tmp_path / "benchmark_tools/results/qfo_corrected_orthomcl_native_sources_20260918.json"
    source_path.parent.mkdir(parents=True)
    source_path.write_text("{}")
    sources = {"status": "corrected_orthomcl_native_sources_prepared_unrun", "threads": 180,
               "pair_workers_planned": 63 if problem == "workers" else 64, "pair_parallel_patch": True,
               "tool_directory": str(tool), "checked_records": [], "originals": [],
               "configured_sources": [module.record(p) for p in tool.iterdir()]}
    if problem == "old_tool_output":
        (tool / "old_run").mkdir()
    elif problem == "old_execution":
        (tool.parent / "inference_execution").mkdir()
    elif problem == "changed_input":
        Path(data["native_inputs"][0]["path"]).write_text("modified")
    monkeypatch.setattr(module, "verify_runtime", lambda root: {})
    monkeypatch.setattr(module, "verify", lambda value: None)
    monkeypatch.setattr(module, "read_frozen", lambda path, sha: data if path == admission_path else sources if path == source_path else {})
    monkeypatch.setattr(module.shutil, "disk_usage", lambda path: SimpleNamespace(free=0 if problem == "space" else 20 * 1024**3))
    def query(argv, **kwargs):
        if argv[0] == "git":
            return "wrong" if problem == "commit" else module.ADMITTER
        return ("JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n"
                "99|COMPLETED|0:0|00:01|bizon|" + ("1" if problem == "allocation" else "2") + "|64G\n")
    monkeypatch.setattr(module.subprocess, "check_output", query)
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: None)
    if problem:
        with pytest.raises((ValueError, FileExistsError)):
            module.preflight(tmp_path, admission_path, "sha", 99)
    else:
        result = module.preflight(tmp_path, admission_path, "sha", 99)
        assert result[0] == data and result[2] == tool
        assert not (tool.parent / "inference_execution").exists()

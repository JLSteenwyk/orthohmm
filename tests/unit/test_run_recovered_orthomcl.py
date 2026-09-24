import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import run_recovered_orthomcl as module


@pytest.mark.parametrize("problem", [None, "allocation", "staging", "exit", "missing_groups",
    "bad_groups", "changed_input", "rewritten_input", "cache", "runtime"])
def test_native_flow(tmp_path, monkeypatch, problem):
    monkeypatch.setenv("SLURM_JOB_ID", "125")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1" if problem == "allocation" else "180")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "921600")
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    output = tmp_path / "native"
    evidence = dict(admission={"fixture": True}, admission_job="124", scheduler={}, accounting="",
                    query_coverage={"failed_queries": 2}, inputs={}, index_validation={})
    checked = []
    monkeypatch.setattr(module, "preflight", lambda *a: (evidence, [], checked, {}, output))
    def sources(tool, inputs, threads):
        assert threads == 180
        tool.mkdir()
        source = tool / "orthomcl.pl"
        source.write_text("fixture")
        return dict(originals=[], configured_sources=[module.record(source)])
    monkeypatch.setattr(module, "write_sources", sources)
    def stage(records, directory, proteins, species, indexes):
        assert (proteins, species) == (984137, 78)
        if problem == "staging":
            raise ValueError("staging failed")
        directory.mkdir()
        staged = {}
        for key, name in dict(bpo="all.bpo", offsets="all_bpo.idx", ranges="all_bpo.se", species="all.gg").items():
            path = directory / name
            path.write_text("fixture")
            staged[key] = module.record(path)
        (directory / "staging.json").write_text("{}")
        return dict(staged=staged)
    monkeypatch.setattr(module, "stage", stage)
    def native(argv, **kwargs):
        assert argv[:3] == ["/usr/bin/time", "-v", "-o"]
        assert kwargs["env"]["ORTHOMCL_PAIR_WORKERS"] == "64"
        assert argv[-6:] == ["--mode", "4", "--bpo_file", str(output / "inputs/all.bpo"),
                            "--gg_file", str(output / "inputs/all.gg")]
        (output / "native.time.txt").write_text("fixture timing")
        if problem != "missing_groups":
            directory = output / "native_tool/run"
            directory.mkdir()
            (directory / "all_orthomcl.out").write_text("fixture groups")
        if problem == "changed_input":
            (output / "inputs/all.bpo").write_text("changed")
        if problem == "rewritten_input":
            path = output / "inputs/all.bpo"
            stat = path.stat()
            module.os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1000000000))
        return SimpleNamespace(returncode=1 if problem == "exit" else 0)
    monkeypatch.setattr(module.subprocess, "run", native)
    monkeypatch.setattr(module, "load_species", lambda *a: {})
    def groups(*args):
        if problem == "bad_groups":
            raise ValueError("bad groups")
        return dict(groups=1, grouped_proteins=2, ungrouped_proteins=984135)
    monkeypatch.setattr(module, "group_counts", groups)
    def caches(*args):
        if problem == "cache":
            raise ValueError("missing cache")
        return []
    monkeypatch.setattr(module, "pair_cache_records", caches)
    def runtime(*args):
        if problem == "runtime":
            raise ValueError("runtime changed")
        return {}
    monkeypatch.setattr(module, "verify_runtime", runtime)
    args = (tmp_path, tmp_path / "admission", "digest", 124, tmp_path / "executor", "commit")
    if problem:
        with pytest.raises(ValueError):
            module.run(*args)
        if problem == "allocation":
            assert not output.exists()
            return
        report = json.loads((output / "report.json").read_text())
        assert report["status"] == "recovered_native_failed"
    else:
        report = module.run(*args)
        assert report["status"] == "recovered_native_exited_zero_pending_admission"
        assert report["content"]["groups"] == 1
    assert report["query_coverage"]["failed_queries"] == 2
    assert report["accuracy_admitted"] is False
    assert report["publication_ready"] is False
    assert report["downstream_execution_authorized"] is False
    assert report["job_id"] == "125"


def test_isolated_help(tmp_path):
    subprocess.run([sys.executable, "-I", "-B", module.__file__, "--help"], cwd=tmp_path,
                   check=True, capture_output=True, text=True)


@pytest.mark.parametrize("problem", [None, "evidence", "runtime", "native_runtime", "changed_source",
                                    "existing", "symlink", "disk"])
def test_preflight(tmp_path, monkeypatch, problem):
    results = tmp_path / "benchmark_tools/results"
    results.mkdir(parents=True)
    source_path = results / "qfo_corrected_orthomcl_native_sources_20260918.json"
    source_path.write_text("{}")
    original = tmp_path / "original.pl"
    original.write_text("original")
    original_record = module.record(original)
    if problem == "changed_source":
        original.write_text("mutated")
    output = tmp_path / "benchmarks/results/qfo_blast_recovery_native_v1"
    output.parent.mkdir(parents=True)
    if problem == "existing":
        output.mkdir()
    elif problem == "symlink":
        output.symlink_to(tmp_path / "absent")
    def runtime(*args):
        if problem == "runtime":
            raise ValueError("wrong Python")
        return {}
    monkeypatch.setattr(module, "verify_runtime", runtime)
    def evidence(*args):
        if problem == "evidence":
            raise ValueError("BPO not admitted")
        return dict(checked_records=[], inputs={"fixture": dict(bytes=1)})
    monkeypatch.setattr(module, "verify_inputs", evidence)
    def frozen(path, digest):
        if path == source_path:
            assert digest == module.SOURCES_SHA
            return dict(originals=[original_record], checked_records=[])
        assert digest in (module.RUNTIME_SHA, module.HELPERS_SHA)
        return {}
    monkeypatch.setattr(module, "read_frozen", frozen)
    def native(*args):
        if problem == "native_runtime":
            raise ValueError("wrong Perl")
    monkeypatch.setattr(module, "verify", native)
    monkeypatch.setattr(module.shutil, "disk_usage", lambda *a: SimpleNamespace(
        free=0 if problem == "disk" else 100 * 1024**3))
    args = (tmp_path, tmp_path / "admission", "digest", 124, tmp_path / "executor", "commit")
    if problem:
        with pytest.raises((ValueError, FileExistsError)):
            module.preflight(*args)
    else:
        result = module.preflight(*args)
        assert result[-1] == output
        assert original_record in result[2]
    if problem not in {"existing", "symlink"}:
        assert not output.exists()

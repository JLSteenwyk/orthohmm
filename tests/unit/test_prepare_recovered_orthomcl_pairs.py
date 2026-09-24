import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools import prepare_recovered_orthomcl_pairs as module


def admission():
    return dict(status="recovered_orthomcl_native_outputs_admitted", conversion_authorized=True,
        admission_job_id="126", node="bizon", allocated_cpus=2, memory_mib=65536,
        accuracy_admitted=False, publication_ready=False, pair_semantics=module.SEMANTICS,
        scheduler=dict(JobIDRaw="125", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="180", ReqMem="900G"),
        content=dict(input_proteins=984137, input_species=78, cross_species_clique_pairs=1,
                     final_groups=1, grouped_proteins=2, ungrouped_input_proteins=984135),
        query_coverage={"failed_queries": 2})


def test_admission():
    value = admission()
    assert module.validate_admission(value, 126) == value["content"]


@pytest.mark.parametrize("key,value", [("status", "corrected_orthomcl_native_outputs_admitted"),
    ("conversion_authorized", False), ("admission_job_id", "125"), ("node", "other"),
    ("allocated_cpus", 1), ("memory_mib", 1024), ("accuracy_admitted", True), ("publication_ready", True),
    ("pair_semantics", "graph_edges")])
def test_contract_rejection(key, value):
    report = admission()
    report[key] = value
    with pytest.raises(ValueError):
        module.validate_admission(report, 126)


@pytest.mark.parametrize("key,value", [("input_proteins", 1), ("input_species", 1),
    ("cross_species_clique_pairs", -1), ("cross_species_clique_pairs", True),
    ("final_groups", -1), ("grouped_proteins", -1), ("ungrouped_input_proteins", 1)])
def test_count_rejection(key, value):
    report = admission()
    report["content"][key] = value
    with pytest.raises(ValueError):
        module.validate_admission(report, 126)


@pytest.mark.parametrize("problem", [None, "audit", "count", "mapping", "changed_input", "runtime_after", "existing"])
def test_flow(tmp_path, monkeypatch, problem):
    root = tmp_path
    report = admission()
    executor = root / "benchmarks/work/executor"
    source = executor / "benchmark_tools/admit_recovered_orthomcl.py"
    source.parent.mkdir(parents=True)
    shutil.copyfile(Path(module.__file__).with_name("admit_recovered_orthomcl.py"), source)
    base = root / "benchmarks/results/qfo_blast_recovery_native_v1"
    native = base / "native_tool/run"
    native.mkdir(parents=True)
    (native / "all_orthomcl.out").write_text("fixture")
    gg = base / "inputs/all.gg"
    gg.parent.mkdir()
    gg.write_text("fixture")
    report.update(source=module.record(source), native_groups=module.record(native / "all_orthomcl.out"),
                  checked_records=[module.record(gg)], outputs=[])
    path = root / "benchmarks/results/qfo_blast_recovery_native_admission_v1/report.json"
    path.parent.mkdir()
    path.write_text(json.dumps(report))
    digest = module.record(path)["sha256"]
    env = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    env.parent.mkdir(parents=True)
    mapping = root / "mapping.json.gz"
    mapping.write_text("fixture")
    env.write_text(json.dumps(dict(reference_files=[module.record(mapping)])))
    monkeypatch.setattr(module, "ENV_SHA", module.record(env)["sha256"])
    monkeypatch.setattr(module, "execution_identity", lambda: {"admission_job_id": "127"})
    monkeypatch.setattr(module, "frozen", lambda *a: None)
    monkeypatch.setattr(module, "completed", lambda job, *a: (report["scheduler"] if job == 125 else {}, ""))
    runtime_calls = []
    def runtime(*args):
        runtime_calls.append(1)
        if problem == "runtime_after" and len(runtime_calls) == 2:
            raise ValueError("runtime changed")
        return {}
    monkeypatch.setattr(module, "verify_runtime", runtime)
    class Owners:
        def __len__(self):
            return 984137

        def values(self):
            return range(78)
    monkeypatch.setattr(module, "load_species", lambda *a: Owners())
    def audit(*args):
        args[-1].write_text("{}")
        content = copy.deepcopy(report["content"])
        if problem == "audit":
            content["final_groups"] += 1
        return dict(content=content, checked_records=[])
    monkeypatch.setattr(module, "audit", audit)
    def write(groups, owners, stream):
        stream.write("A\tB\n")
        return dict(total_pairs=2 if problem == "count" else 1, final_groups=1,
                    grouped_proteins=2, ungrouped_input_proteins=984135)
    monkeypatch.setattr(module, "write_pairs", write)
    monkeypatch.setattr(module, "load_mapping", lambda *a: {})
    def filter_pairs(source, output, mapping):
        output.write_text(source.read_text())
        if problem == "changed_input":
            gg.write_text("mutated")
        return 1, 0 if problem == "mapping" else 1
    monkeypatch.setattr(module, "filter_pairs", filter_pairs)
    destination = root / "benchmarks/results/qfo_blast_recovery_pairs_v1"
    if problem == "existing":
        destination.mkdir()
    if problem:
        with pytest.raises((ValueError, FileExistsError)):
            module.prepare(root, path, digest, 126, executor, "commit")
        if problem == "existing":
            assert not (destination / "results.json").exists()
            return
        result = json.loads((destination / "results.json").read_text())
        assert result["status"] == "recovered_pair_conversion_failed"
        assert not (destination / "pairs.tsv").exists()
    else:
        result = module.prepare(root, path, digest, 126, executor, "commit")
        assert result["status"] == "recovered_orthomcl_pairs_prepared_unscored"
        assert result["total_pairs"] == result["retained_pairs"] == 1
        assert (destination / "pairs.qfo.tsv").read_text() == "A\tB\n"
    assert result["accuracy_evaluated"] is False
    assert result["publication_ready"] is False
    assert result["query_coverage"]["failed_queries"] == 2


def test_isolated_cli_and_pins(tmp_path):
    for name, digest in (("admit_recovered_orthomcl.py", module.ADMITTER_SHA),
                         ("prepare_qfo_corrected_orthomcl_pairs.py", module.CONVERTER_SHA)):
        assert module.record(Path(module.__file__).with_name(name))["sha256"] == digest
    subprocess.run([sys.executable, "-I", "-B", module.__file__, "--help"], cwd=tmp_path,
                   check=True, capture_output=True, text=True)

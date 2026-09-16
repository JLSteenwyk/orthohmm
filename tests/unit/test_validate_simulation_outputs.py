import json

import pytest

from benchmark_tools.benchmark_production import file_record
from benchmark_tools.validate_simulation_outputs import validate_graph_weights, validate_orthofinder, validate_orthohmm, verify_process


def inventory(root, extras=()):
    return [dict(file_record(p, p.parent), absolute_path=str(p))
            for p in [*sorted(p for p in root.rglob("*") if p.is_file()), *extras]]


def test_process_requires_command_success_and_exact_output_set(tmp_path):
    root = tmp_path / "out"
    root.mkdir()
    (root / "result").write_text("result")
    config = {"output": str(root), "argv": ["tool"]}
    record = {"argv": ["tool"], "status": "process_succeeded", "exit_code": 0, "outputs": inventory(root)}
    assert verify_process(config, record) == {root / "result"}
    record["exit_code"] = 1
    with pytest.raises(ValueError, match="did not succeed"):
        verify_process(config, record)
    record["exit_code"] = 0
    (root / "extra").write_text("extra")
    with pytest.raises(ValueError, match="file set changed"):
        verify_process(config, record)


def test_orthofinder_native_completion_and_version(tmp_path):
    root = tmp_path / "of"
    inputs = root / "input"
    inputs.mkdir(parents=True)
    fasta = inputs / "a.fasta"
    fasta.write_text(">a\nAAA\n")
    input_records = [dict(file_record(fasta, inputs), absolute_path=str(fasta))]
    log = root / "Log.txt"
    log.write_text("2026 : Started OrthoFinder version 3.1.5\nCommand Line: tool -f input\n2026 : OrthoFinder run completed\n")
    (root / "OrthoFinder_graph.txt").write_text("0 1:0.9 $\n")
    config = {"argv": ["tool", "-f", "input"], "output": str(root), "copy_inputs_to": str(inputs)}
    record = {"argv": config["argv"], "status": "process_succeeded", "exit_code": 0, "outputs": inventory(root)}
    assert "verified" in validate_orthofinder(config, record, input_records)["completion"]
    log.write_text(log.read_text().replace("3.1.5", "2.5.5"))
    record["outputs"] = inventory(root)
    with pytest.raises(ValueError, match="completion not confirmed"):
        validate_orthofinder(config, record, input_records)
    log.write_text("2026 : Started OrthoFinder version 3.1.5\nCommand Line: tool -f input\n")
    record["outputs"] = inventory(root)
    with pytest.raises(ValueError, match="completion not confirmed"):
        validate_orthofinder(config, record, input_records)


@pytest.mark.parametrize("weight", ["nan", "inf", "-inf"])
def test_nonfinite_clustering_graph_rejected(tmp_path, weight):
    graph = tmp_path / "graph"
    graph.write_text(f"(mclmatrix\nbegin\n0 1:{weight} $\n)\n")
    with pytest.raises(ValueError, match="Nonfinite"):
        validate_graph_weights(graph)


def test_orthohmm_requires_native_completion_and_frozen_sources(tmp_path):
    root = tmp_path / "out"
    root.mkdir()
    native = root / "orthohmm_orthogroups.txt"
    native.write_text("OG1: a b\n")
    source = tmp_path / "source"
    source.mkdir()
    module = source / "core.py"
    module.write_text("pass\n")
    fasta = tmp_path / "a.fasta"
    fasta.write_text(">a\nAAA\n")
    inputs = [dict(file_record(fasta, tmp_path), absolute_path=str(fasta))]
    manifest = {"core_root": str(source), "core_commit": "frozen",
                "core_sources": [dict(file_record(module, source), absolute_path=str(module))]}
    metrics = tmp_path / "metrics.json"
    data = {"status": "complete", "harness": {"exit_code": 0, "git_commit": "frozen",
            "input_manifest": [file_record(fasta, tmp_path)], "source_manifest": [file_record(module, source)],
            "output_manifest": [file_record(native, root)]}}
    metrics.write_text(json.dumps(data))
    config = {"output": str(root), "metrics": str(metrics), "argv": ["tool"]}
    record = {"status": "process_succeeded", "exit_code": 0, "argv": ["tool"], "outputs": inventory(root, [metrics])}
    assert validate_orthohmm("orthohmm_high_sensitivity", config, record, inputs, manifest)["native"] == str(native)
    data["status"] = "failed"
    metrics.write_text(json.dumps(data))
    record["outputs"] = inventory(root, [metrics])
    with pytest.raises(ValueError, match="completion not confirmed"):
        validate_orthohmm("orthohmm_high_sensitivity", config, record, inputs, manifest)

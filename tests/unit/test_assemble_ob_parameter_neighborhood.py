import copy
import json

import pytest

from benchmark_tools import assemble_ob_parameter_neighborhood as assembler


def test_pending_native_run_prevents_reference_access(tmp_path, monkeypatch):
    def fail(*args):
        raise ValueError("Pending native run")
    monkeypatch.setattr(assembler, "admit_panel", fail)
    monkeypatch.setattr(assembler, "load_reference_snapshot", lambda *a: pytest.fail("Premature scoring"))
    with pytest.raises(ValueError, match="Pending"):
        assembler.assemble(tmp_path, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_no_overwrite(tmp_path, monkeypatch):
    monkeypatch.setattr(assembler, "admit_panel", lambda *a: pytest.fail("Overwriting"))
    with pytest.raises(FileExistsError):
        assembler.assemble(tmp_path, tmp_path)


def test_admission_uses_each_frozen_panel(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(assembler, "validate_baseline", lambda root, index: (root, index))
    def variant(root, index, cpm):
        calls.append((index, cpm))
        return (root, index, cpm)
    monkeypatch.setattr(assembler, "validate_variant", variant)
    panel = assembler.admit_panel(tmp_path)
    assert panel["control"] == (tmp_path, 3)
    assert calls == [(0, True), (1, True), (0, False), (1, False), (2, False), (3, False)]
    assert tuple(panel["variants"]) == assembler.VARIANTS


@pytest.mark.parametrize("problem", [None, "official", "mutation", "coverage", "control", "missing", "identity", "native"])
def test_full_panel_scoring(tmp_path, monkeypatch, problem):
    fasta = tmp_path / "species.fa"
    fasta.write_text(">a\nACD\n>b\nACD\n>c\nACD\n")
    prediction = tmp_path / "groups.tsv"
    prediction.write_text("root_hog\tsource_family\tgenes\nRootHOG0\tFamily0\ta,b\n" +
                          ("" if problem == "coverage" else "RootHOG1\tFamily1\tc\n"))
    metrics = tmp_path / "metrics.json"
    metrics.write_text(json.dumps({"outputs": {"root_hogs": assembler.file_provenance(prediction)},
                                  "rss_measurement": "sampled_sum_of_linux_proc_tree_rss",
                                  "wall_s": 10, "user_cpu_s": 20, "system_cpu_s": 1,
                                  "peak_process_tree_rss_bytes": 1024}))
    admission = {"native_metrics": assembler.file_provenance(metrics), "status": "native_group_output_verified",
                 "native_outputs_validated": True, "accuracy_evaluated": False}
    panel = {"control": copy.deepcopy(admission), "variants": {
        label: {"status": "candidate_variant_native_validated_unscored", "variant": label,
                "accuracy_evaluated": False, "native_validation": copy.deepcopy(admission)}
        for label in assembler.VARIANTS}}
    first = panel["variants"][assembler.VARIANTS[0]]
    if problem == "control":
        panel["control"]["native_outputs_validated"] = False
    elif problem == "missing":
        panel["variants"].pop(assembler.VARIANTS[0])
    elif problem == "identity":
        first["variant"] = "another_arm"
    elif problem == "native":
        first["native_validation"]["native_outputs_validated"] = False
    monkeypatch.setattr(assembler, "admit_panel", lambda *a: panel)
    monkeypatch.setattr(assembler, "read_frozen", lambda *a: {"fasta_inputs": [assembler.file_provenance(fasta)]})
    official = tmp_path / "official.py"
    official.write_text("# Fixture\n")
    def references(*args):
        assert problem not in {"control", "missing", "identity", "native", "coverage"}
        return {"RefOG001.txt": {"a", "b"}}, {}, official, []
    monkeypatch.setattr(assembler, "load_reference_snapshot", references)
    def official_score(*args):
        if problem == "mutation":
            fasta.write_text(">a\nAAA\n>b\nACD\n>c\nACD\n")
        return {"f_score": 99 if problem == "official" else 100,
                "precision": 100, "recall": 100, "exact_refogs": 1}
    monkeypatch.setattr(assembler, "run_official_benchmark", official_score)
    output = tmp_path / "out"
    if problem:
        with pytest.raises(ValueError):
            assembler.assemble(tmp_path, output)
        assert not (output / "results.json").exists()
        return
    result = assembler.assemble(tmp_path, output)
    assert result["replicates"] == 20000 and result["seed"] == 20260918
    assert result["planned_endpoints"] == 18
    assert len(result["comparisons"]) == 6
    assert result["failed_variants"] == {} and result["publication_ready"] is False
    for label in assembler.VARIANTS:
        for metric in result["comparisons"][label]["metrics"].values():
            assert metric["difference_percentage_points"] == 0
            assert metric["bonferroni_percentile_ci"] == [0, 0]
        assert result["coverage_resources"][label]["assigned_genes"] == 3
    assert "not end-to-end" in (output / "results.md").read_text()

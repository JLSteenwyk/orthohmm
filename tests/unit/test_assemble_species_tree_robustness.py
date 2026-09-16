import copy
import json

import pytest

from benchmark_tools import assemble_species_tree_robustness as assembler


def test_native_failure_prevents_reference_access(tmp_path, monkeypatch):
    def fail(*args):
        raise ValueError("Incomplete panel")
    monkeypatch.setattr(assembler, "validate", fail)
    monkeypatch.setattr(assembler, "load_reference_snapshot", lambda *a: pytest.fail("Read labels before admission"))
    with pytest.raises(ValueError, match="Incomplete"):
        assembler.assemble(tmp_path, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_no_overwrite_before_admission(tmp_path, monkeypatch):
    monkeypatch.setattr(assembler, "validate", lambda *a: pytest.fail("Overwriting output"))
    with pytest.raises(FileExistsError):
        assembler.assemble(tmp_path, tmp_path)


@pytest.mark.parametrize("change", ["missing", "extra"])
def test_exact_comparison_inventory(change):
    scores = {name: {} for name in (assembler.BASELINE, *assembler.PERTURBATIONS)}
    if change == "missing":
        scores.pop(assembler.PERTURBATIONS[0])
    else:
        scores["extra"] = {}
    with pytest.raises(ValueError, match="exactly"):
        assembler.bootstrap(scores)


@pytest.mark.parametrize("problem", [None, "official", "mutation", "coverage", "control", "missing", "native"])
def test_panel_scoring_admission_and_fixed_endpoints(tmp_path, monkeypatch, problem):
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
    panel = {"status": "native_panel_validated_unscored", "accuracy_evaluated": False,
             "control_validation": {"status": "equivalent", "native_validation": admission},
             "variants": {label: {"native_validation": copy.deepcopy(admission), "rooted_rf_clade_distance": 2 if i < 3 else 4}
                          for i, label in enumerate(assembler.PERTURBATIONS)}}
    if problem == "control":
        panel["control_validation"]["status"] = "not_equivalent"
    elif problem == "missing":
        panel["variants"].pop(assembler.PERTURBATIONS[0])
    elif problem == "native":
        panel["variants"][assembler.PERTURBATIONS[0]]["native_validation"]["native_outputs_validated"] = False
    monkeypatch.setattr(assembler, "validate", lambda *a: panel)
    monkeypatch.setattr(assembler, "read_frozen", lambda *a: {"fasta_inputs": [assembler.file_provenance(fasta)]})
    official = tmp_path / "official.py"
    official.write_text("# Fixture\n")
    def references(*args):
        assert problem not in {"control", "missing", "native", "coverage"}
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
    assert "18 reported" in result["multiplicity"]
    assert len(result["comparisons"]) == 6
    for label in assembler.PERTURBATIONS:
        contrast = result["comparisons"][label]
        assert contrast["family_f1_ties"] == 1
        for metric in contrast["metrics"].values():
            assert metric["difference_percentage_points"] == 0
            assert metric["bonferroni_percentile_ci"] == [0, 0]
        assert result["coverage_resources"][label]["assigned_genes"] == 3
    assert json.loads((output / "results.json").read_text())["publication_ready"] is False
    assert "not end-to-end" in (output / "results.md").read_text()

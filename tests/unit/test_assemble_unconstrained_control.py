import json

import pytest

from benchmark_tools import assemble_unconstrained_control as assembler


def test_native_failure_prevents_reference_access(tmp_path, monkeypatch):
    def fail(*args):
        raise ValueError("Still running or invalid native output")

    def forbidden(*args):
        pytest.fail("Read labels before native admission")

    monkeypatch.setattr(assembler, "validate_unconstrained", fail)
    monkeypatch.setattr(assembler, "load_reference_snapshot", forbidden)
    monkeypatch.setattr(assembler, "read_frozen", forbidden)
    output = tmp_path / "result"
    with pytest.raises(ValueError, match="Still running"):
        assembler.assemble(tmp_path, output)
    assert not output.exists()


def test_no_overwrite_even_before_validation(tmp_path, monkeypatch):
    monkeypatch.setattr(assembler, "validate_unconstrained", lambda *a: pytest.fail("Overwriting output"))
    with pytest.raises(FileExistsError):
        assembler.assemble(tmp_path, tmp_path)


def test_rejects_extra_comparison():
    with pytest.raises(ValueError, match="exactly"):
        assembler.bootstrap({assembler.BASELINE: {}, assembler.DIAGNOSTIC: {}, "extra": {}})


@pytest.mark.parametrize("problem", [None, "official", "mutation", "coverage"])
def test_pair_assembly_crosschecks_and_fixed_bootstrap(tmp_path, monkeypatch, problem):
    fasta = tmp_path / "species.fa"
    fasta.write_text(">a\nACD\n>b\nACD\n>c\nACD\n")
    metrics = tmp_path / "metrics.json"
    metrics.write_text(json.dumps({"rss_measurement": "sampled_sum_of_linux_proc_tree_rss",
                                  "wall_s": 10, "user_cpu_s": 20, "system_cpu_s": 1,
                                  "peak_process_tree_rss_bytes": 1024}))
    native = {"native_metrics": assembler.file_provenance(metrics)}
    admitted = []

    def validate(*args):
        admitted.append(args)
        return native

    monkeypatch.setattr(assembler, "validate_unconstrained", validate)
    monkeypatch.setattr(assembler, "validate_constrained", validate)
    prepared = {"fasta_inputs": [assembler.file_provenance(fasta)]}
    monkeypatch.setattr(assembler, "read_frozen", lambda *a: prepared)
    cells = []
    for label in (assembler.BASELINE, assembler.DIAGNOSTIC):
        path = tmp_path / f"{label}.tsv"
        path.write_text("root_hog\tsource_family\tgenes\nRootHOG0\tFamily0\ta,b\n" +
                        ("" if problem == "coverage" else "RootHOG1\tFamily1\tc\n"))
        cells.append({"label": label, "prediction": str(path)})
    monkeypatch.setattr(assembler, "select_cell", lambda *a: (cells[0], tmp_path, tmp_path))
    monkeypatch.setattr(assembler, "unconstrained_cell", lambda *a: cells[1])
    official = tmp_path / "official.py"
    official.write_text("# Test fixture\n")

    def references(*args):
        assert len(admitted) == 2
        return {"RefOG001.txt": {"a", "b"}}, {}, official, []

    monkeypatch.setattr(assembler, "load_reference_snapshot", references)

    def official_score(*args):
        if problem == "mutation":
            fasta.write_text(">a\nAAA\n>b\nACD\n>c\nACD\n")
        return {"f_score": 99 if problem == "official" else 100,
                "precision": 100, "recall": 100, "exact_refogs": 1}

    monkeypatch.setattr(assembler, "run_official_benchmark", official_score)
    output = tmp_path / "assembled"
    if problem:
        with pytest.raises(ValueError):
            assembler.assemble(tmp_path, output)
        assert not (output / "results.json").exists()
        return
    result = assembler.assemble(tmp_path, output)
    assert result["replicates"] == 20000
    assert result["seed"] == 20260918
    assert "3 reported" in result["multiplicity"]
    assert len(result["comparisons"]) == 1
    assert result["comparisons"][assembler.DIAGNOSTIC]["family_f1_ties"] == 1
    for metric in result["comparisons"][assembler.DIAGNOSTIC]["metrics"].values():
        assert metric["difference_percentage_points"] == 0
        assert metric["bonferroni_percentile_ci"] == [0, 0]
    assert result["coverage_resources"][assembler.DIAGNOSTIC]["assigned_genes"] == 3
    assert json.loads((output / "results.json").read_text())["publication_ready"] is False
    assert "not a ninth factorial cell" in (output / "results.md").read_text()

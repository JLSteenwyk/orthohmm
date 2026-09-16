from pathlib import Path
import json

import pytest

from benchmark_tools import assemble_orthobench_factorial as assembler


def accounting(state="FAILED"):
    return "JobID|State|ExitCode\n" + "".join(f"21248_{i}|{state}|1:0\n" for i in range(4))


def test_accepts_complete_terminal_task_set_without_relabeling():
    rows = assembler.require_terminal_array(accounting())
    assert len(rows) == 4
    assert all(r["State"] == "FAILED" for r in rows)


@pytest.mark.parametrize("state", ["RUNNING", "PENDING", "COMPLETING", "UNKNOWN"])
def test_rejects_nonterminal_tasks(state):
    with pytest.raises(ValueError):
        assembler.require_terminal_array(accounting(state))


@pytest.mark.parametrize("change", ["duplicate", "missing"])
def test_rejects_ambiguous_task_set(change):
    text = accounting()
    text = text + "21248_0|FAILED|1:0\n" if change == "duplicate" else text.replace("21248_3|FAILED|1:0\n", "")
    with pytest.raises(ValueError):
        assembler.require_terminal_array(text)


def test_running_array_never_reaches_validation_or_scoring(tmp_path, monkeypatch):
    monkeypatch.setattr(assembler.subprocess, "check_output", lambda *a, **k: accounting("RUNNING"))
    def forbidden(*args, **kwargs):
        pytest.fail("Read inference or reference outcomes before terminal gate")
    monkeypatch.setattr(assembler, "validate", forbidden)
    monkeypatch.setattr(assembler, "read_frozen", forbidden)
    monkeypatch.setattr(assembler, "score_partition", forbidden)
    output = tmp_path / "results"
    with pytest.raises(ValueError, match="no partial scoring"):
        assembler.assemble(Path.cwd(), output)
    assert not output.exists()


def test_official_rounding_is_percentage_not_fraction():
    exact = {"f_score": 74.106, "precision": 81.770, "recall": 67.755, "exact_refogs": 20}
    printed = {"f_score": 74.1, "precision": 81.8, "recall": 67.8, "exact_refogs": 20}
    assembler.compare_official(exact, printed)


@pytest.mark.parametrize("change", ["f_score", "precision", "recall", "exact_refogs", "missing", "nan"])
def test_rejects_scorer_disagreement(change):
    exact = {"f_score": 70.0, "precision": 70.0, "recall": 70.0, "exact_refogs": 20}
    printed = exact.copy()
    if change == "missing":
        del printed["precision"]
    elif change == "nan":
        printed["recall"] = float("nan")
    else:
        printed[change] += 1
    with pytest.raises(ValueError):
        assembler.compare_official(exact, printed)


def test_coverage_distinguishes_singletons_paralogs_and_multispecies_groups():
    groups = [{"a", "b"}, {"c", "d"}, {"e"}]
    species = {"a": "s1", "b": "s1", "c": "s1", "d": "s2", "e": "s2"}
    result = assembler.coverage_and_resources(groups, species)
    assert result["assigned_genes"] == 5
    assert result["singleton_groups"] == 1
    assert result["genes_in_nonsingleton_groups"] == 4
    assert result["multispecies_groups"] == 1
    assert result["genes_in_multispecies_groups"] == 2
    assert result["wall_s"] is None


def resource_metrics():
    return {"rss_measurement": "sampled_sum_of_linux_proc_tree_rss", "wall_s": 10,
            "user_cpu_s": 20, "system_cpu_s": 5, "peak_process_tree_rss_bytes": 2**30}


def test_resource_accounting_uses_cpu_seconds_and_keeps_bytes():
    result = assembler.coverage_and_resources([{"a"}], {"a": "s1"}, resource_metrics())
    assert result["mean_cpu_cores"] == 2.5
    assert result["peak_process_tree_rss_bytes"] == 2**30
    assert "shared node" in result["resource_scope"]


@pytest.mark.parametrize("key,value", [("wall_s", 0), ("user_cpu_s", -1), ("system_cpu_s", float("nan")),
                                       ("peak_process_tree_rss_bytes", True), ("rss_measurement", "maxrss")])
def test_rejects_invalid_or_incommensurate_resources(key, value):
    metrics = resource_metrics()
    metrics[key] = value
    with pytest.raises(ValueError):
        assembler.coverage_and_resources([{"a"}], {"a": "s1"}, metrics)


def test_missing_input_gene_is_not_reported_as_full_coverage():
    with pytest.raises(ValueError):
        assembler.coverage_and_resources([{"a"}], {"a": "s1", "b": "s2"})


def test_execution_table_labels_unmeasured_resources():
    row = assembler.coverage_and_resources([{"a"}], {"a": "s1"})
    report = assembler.render_execution({"coverage_resources": {label: row for label in assembler.CELLS}})
    assert report.count("NA | NA | NA") == 8
    assert "not zero cost" in report
    assert "not orthology accuracy" in report


@pytest.mark.parametrize("disagree", [False, True])
def test_eight_cell_assembly_writes_reports_only_after_crosschecks(tmp_path, monkeypatch, disagree):
    fasta = tmp_path / "species.fa"
    fasta.write_text(">a\nACD\n>b\nACD\n>c\nACD\n")
    metrics = tmp_path / "metrics.json"
    metrics.write_text(json.dumps(resource_metrics()))
    ref = tmp_path / "RefOG001.txt"
    ref.write_text("a\nb\n")
    official = tmp_path / "official.py"
    official.write_text("# Golden expected values supplied by this unit test\n")
    cells, arms = [], {}
    for label in assembler.CELLS:
        p, c, r = (bool(int(part[1])) for part in label.split("_"))
        path = tmp_path / (label + ".txt")
        path.write_text("root_hog\tsource_family\tgenes\nRootHOG0000000\tFamily0000000\ta\n"
                        "RootHOG0000001\tFamily0000000\tb\nRootHOG0000002\tFamily0000001\tc\n"
                        if r else "a b\nc\n")
        cells.append({"label": label, "profile_expansion": p, "candidate_expansion": c,
                      "reconciliation": r, "prediction": str(path)})
        if not r:
            arms[f"p{int(p)}_c{int(c)}"] = {"candidate_partition": assembler.file_provenance(path)}
    prepared = {"cells": cells, "candidate_arms": arms, "fasta_inputs": [{"path": str(fasta)}]}
    monkeypatch.setattr(assembler.subprocess, "check_output", lambda *a, **k: accounting())
    monkeypatch.setattr(assembler, "read_frozen", lambda *a: prepared)
    validated = []

    def validate(*args):
        validated.append(args[1])
        return {"native_metrics": assembler.file_provenance(metrics)}

    monkeypatch.setattr(assembler, "validate", validate)
    monkeypatch.setattr(assembler, "select_cell", lambda *a: (cells[1], tmp_path, tmp_path))
    monkeypatch.setattr(assembler, "verify_prepared", lambda *a: None)
    monkeypatch.setattr(assembler, "load_reference_snapshot", lambda *a: (
        {"RefOG001.txt": {"a", "b"}}, {}, official, [assembler.file_provenance(ref)]))
    crosschecked = []

    def golden_official(script, path):
        assert validated == [0, 1, 2, 3]
        lines = path.read_text().splitlines()
        reconciled = path.stem.endswith("r1")
        assert lines == (["a", "b", "c"] if reconciled else ["a b", "c"])
        crosschecked.append(path.stem)
        score = 0.0 if reconciled else 100.0
        return {"f_score": score - int(disagree), "precision": score, "recall": score,
                "exact_refogs": 0 if reconciled else 1}

    monkeypatch.setattr(assembler, "run_official_benchmark", golden_official)
    output = tmp_path / "assembled"
    if disagree:
        with pytest.raises(ValueError, match="Official scorer disagreement"):
            assembler.assemble(tmp_path, output)
        assert not (output / "results.json").exists()
        assert not (output / "results.md").exists()
        return
    result = assembler.assemble(tmp_path, output)
    assert crosschecked == list(assembler.CELLS)
    assert len(result["comparisons"]) == 12
    assert result["draws_generated"] == 20000
    assert result["scores"]["p0_c0_r0"]["f_score"] == 100.0
    assert result["scores"]["p0_c0_r1"]["f_score"] == 0.0
    assert result["coverage_resources"]["p0_c0_r0"]["wall_s"] is None
    assert json.loads((output / "results.json").read_text())["publication_ready"] is False
    assert "Coverage And Incremental Resources" in (output / "results.md").read_text()
    with pytest.raises(FileExistsError):
        assembler.assemble(tmp_path, output)

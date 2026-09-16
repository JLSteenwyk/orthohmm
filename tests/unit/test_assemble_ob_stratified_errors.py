import copy
import json

import pytest

from benchmark_tools import assemble_ob_stratified_errors as assembler
from benchmark_tools.prepare_ob_error_strata import assign_strata


def inputs():
    records = [{"refog": f"r{i}", "genes": 2+i, "true_positive": (2+i)*(1+i)//2,
                "false_positive": 0, "false_negative": 0} for i in range(6)]
    scores = {name: {"refog_records": copy.deepcopy(records)} for name in assembler.METHODS}
    categories = {d + ":" + c: [] for d, cats in assembler.CATEGORIES.items() for c in cats}
    for dimension, cats in assembler.CATEGORIES.items():
        categories[dimension + ":" + cats[0]] = [row["refog"] for row in records]
    categories["composition:concentrated"] = ["r0"]
    categories["composition:not_concentrated"] = [f"r{i}" for i in range(1, 6)]
    return scores, categories


def test_all_endpoints_retained_with_explicit_small_and_empty_bins():
    scores, categories = inputs()
    result = assembler.stratified_statistics(scores, categories)
    assert len(result) == 14
    assert sum(len(c["metrics"]) for row in result.values() for c in row["comparisons"].values()) == 84
    large = result["composition:not_concentrated"]
    assert large["status"] == "paired_bootstrap"
    assert large["bootstrap"]["replicates"] == 20000
    assert large["bootstrap"]["seed"] == 20260918
    assert "84 planned" in large["bootstrap"]["multiplicity"]
    small = result["composition:concentrated"]
    assert small["status"] == "descriptive_only_lt5"
    assert small["point_estimates_percent"][assembler.BASELINE]["f_score"] == 100
    empty = result["composition:missing"]
    assert empty["status"] == "empty_nonestimable"
    assert empty["point_estimates_percent"][assembler.BASELINE] is None
    for row in (small, empty):
        assert row["multiplicity_endpoints"] == 84
        for comparison in row["comparisons"].values():
            for value in comparison["metrics"].values():
                assert value["paired_percentile_ci"] is None
                assert value["bonferroni_percentile_ci"] is None


@pytest.mark.parametrize("problem", ["method", "bin", "duplicate_family", "family_size"])
def test_incomplete_or_incompatible_panels_rejected(problem):
    scores, categories = inputs()
    if problem == "method":
        scores.pop(assembler.COMPARATORS[0])
    elif problem == "bin":
        categories.pop("composition:missing")
    elif problem == "duplicate_family":
        categories["composition:missing"] = ["r0"]
    elif problem == "family_size":
        scores[assembler.COMPARATORS[0]]["refog_records"][0].update(genes=3, true_positive=3)
    with pytest.raises(ValueError):
        assembler.stratified_statistics(scores, categories)


def test_failed_feature_admission_prevents_outcome_access(tmp_path, monkeypatch):
    monkeypatch.setattr(assembler, "read_frozen", lambda path, sha: {} if path.name.startswith("ob_error_strata") else pytest.fail("Read outcomes"))
    monkeypatch.setattr(assembler, "validate_strata", lambda *a: None)
    def fail(*args):
        raise ValueError("Incomplete native admission")
    monkeypatch.setattr(assembler, "recheck_strata", fail)
    with pytest.raises(ValueError, match="Incomplete"):
        assembler.assemble(tmp_path, tmp_path / "out")
    assert not (tmp_path / "out/results.json").exists()


@pytest.mark.parametrize("problem", [None, "official", "mutation", "statistics", "unknown_gene"])
def test_native_parsing_and_full_reference_crosschecks(tmp_path, monkeypatch, problem):
    references = {f"RefOG{i:03d}.txt": {f"g{i}a", f"g{i}b"} for i in range(1, 71)}
    families = {name: {"genes": 2, "residue_lengths": [100, 100], "composition_flags": [False, False],
                      "species_counts": {"s": 2}, "mean_pairwise_identity": .5} for name in references}
    frozen = {**assign_strata(families), "families": families, "status": "family_strata_prepared_unscored", "accuracy_evaluated": False}
    fasta = tmp_path / "species.fa"
    fasta.write_text("".join(f">{gene}\nACD\n" for group in references.values() for gene in sorted(group)))
    prepared = {"fasta_inputs": [assembler.file_provenance(fasta)]}
    groups = list(references.values())
    sources = {}
    for method in assembler.METHODS:
        path = tmp_path / (method + ".txt")
        if method == assembler.BASELINE:
            content = "".join(f"OG{i}: {' '.join(sorted(group))}\n" for i, group in enumerate(groups))
        elif method == assembler.COMPARATORS[1]:
            content = "root_hog\tsource_family\tgenes\n" + "".join(
                f"H{i}\tF{i}\t{','.join(sorted(group))}\n" for i, group in enumerate(groups))
        else:
            content = "".join(" ".join(sorted(group)) + "\n" for group in groups)
            if problem == "unknown_gene":
                content += "unknown_gene\n"
        path.write_text(content)
        sources[method] = assembler.file_provenance(path)
    score = assembler.score_partition(groups, references, {})
    snapshot = {"scores": {method: copy.deepcopy(score) for method in assembler.METHODS}, "inputs": {"predictions": sources}}
    if problem == "statistics":
        snapshot["scores"][assembler.BASELINE]["refog_records"][0].update(true_positive=0., false_negative=1.)
    def read(path, sha):
        if path.name.startswith("ob_error_strata"):
            return frozen
        if path.name.startswith("orthobench_paired"):
            return snapshot
        return prepared
    monkeypatch.setattr(assembler, "read_frozen", read)
    def recheck(root, output):
        output.mkdir(parents=True)
        return frozen
    monkeypatch.setattr(assembler, "recheck_strata", recheck)
    official = tmp_path / "official.py"
    official.write_text("# fixture\n")
    monkeypatch.setattr(assembler, "load_reference_snapshot", lambda *a: (references, {}, official, []))
    real_provenance = assembler.file_provenance
    monkeypatch.setattr(assembler, "file_provenance", lambda path: real_provenance(path) if path.exists() else {"path": str(path)})
    def official_score(*args):
        if problem == "mutation":
            fasta.write_text(fasta.read_text().replace("ACD", "AAA"))
        return {"f_score": 99 if problem == "official" else 100, "precision": 100, "recall": 100, "exact_refogs": 70}
    monkeypatch.setattr(assembler, "run_official_benchmark", official_score)
    output = tmp_path / "out"
    if problem:
        with pytest.raises(ValueError):
            assembler.assemble(tmp_path, output)
        assert not (output / "results.json").exists()
    else:
        result = assembler.assemble(tmp_path, output)
        assert result["status"] == "stratified_analysis_complete"
        assert result["multiplicity_endpoints"] == 84
        assert all(row["observed_prediction_genes"] == 140 for row in result["native_coverage"].values())
        assert json.loads((output / "results.json").read_text())["publication_ready"] is False
        assert "Empty bins are not scored as zero" in (output / "results.md").read_text()

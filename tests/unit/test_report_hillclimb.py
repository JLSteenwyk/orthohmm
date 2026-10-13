import json
from pathlib import Path

import pytest

from benchmark_tools import report_hillclimb as reporting
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("value", [-1, 1.1, float("nan"), float("inf"), True, "0.8", None])
def test_invalid_score_values_are_rejected(value):
    with pytest.raises(ValueError, match="Invalid"):
        reporting.endpoint_vector("OrthoBench", {"OrthoBench F1": value})


@pytest.mark.parametrize("scores", [{}, {"unknown": .8}, {"OrthoBench F1": .8, "extra": .8}])
def test_endpoint_identity_is_strict(scores):
    with pytest.raises(ValueError, match="endpoint"):
        reporting.endpoint_vector("OrthoBench", scores)


def test_zero_is_measured_not_missing():
    assert reporting.endpoint_vector("OrthoBench", {"OrthoBench F1": 0}) == {ENDPOINTS[0]: 0}


def test_complete_vector_screen_requires_primary_gain_and_no_large_decline():
    baseline = dict.fromkeys(ENDPOINTS, .8)
    primary = {**baseline, ENDPOINTS[0]: .81}
    secondary = {**baseline, ENDPOINTS[4]: .81}
    tradeoff = {**primary, ENDPOINTS[5]: .79}
    partial = {k: v for k, v in primary.items() if k != ENDPOINTS[-1]}
    result = reporting.screen({"control": baseline, "primary": primary, "secondary": secondary,
                               "tradeoff": tradeoff, "partial": partial}, "control")
    assert result["comparisons"]["primary"]["passes_point_estimate_screen"]
    assert not result["comparisons"]["secondary"]["passes_point_estimate_screen"]
    assert not result["comparisons"]["tradeoff"]["passes_point_estimate_screen"]
    assert not result["comparisons"]["partial"]["passes_point_estimate_screen"]
    assert result["pareto"]["excluded_incomplete"]["partial"] == [ENDPOINTS[-1]]


def test_incomplete_matched_control_cannot_qualify_complete_variant():
    scores = dict.fromkeys(ENDPOINTS, .8)
    result = reporting.screen({"control": {ENDPOINTS[0]: .7}, "variant": scores}, "control")
    assert not result["comparisons"]["variant"]["passes_point_estimate_screen"]


def test_retained_schema_and_duplicate_cells_are_rejected():
    row = {"dataset": "OrthoBench", "key": "tool", "label": "Tool", "scores": {"weighted_refog_F1": .8},
           "declared_version": "1", "output_semantics": "groups", "resources": [], "gaps": []}
    result = reporting.retained_methods({"rows": [row]})
    assert result["tool"]["scores"] == {ENDPOINTS[0]: .8}
    with pytest.raises(ValueError, match="Duplicate"):
        reporting.retained_methods({"rows": [row, row]})
    with pytest.raises(ValueError, match="schema"):
        reporting.retained_methods({"rows": [{**row, "scores": {"F1": .8}}]})


def cell(score=.8, identity="same"):
    return {"method_key": "membership:control_direct", "dataset": "ThreeKingdoms",
            "status": "independent_content_readback_verified", "scores": {ENDPOINTS[-1]: score},
            "products": {name: {"sha256": identity} for name in reporting.PRODUCTS}}


def test_duplicate_control_requires_matching_scores_and_native_outputs():
    result = reporting.combine_cells({}, [cell(), cell()])
    assert result["membership:control_direct"]["scores"] == {ENDPOINTS[-1]: .8}
    assert len(result["membership:control_direct"]["cells"]) == 2
    with pytest.raises(ValueError, match="Conflicting"):
        reporting.combine_cells({}, [cell(), cell(.81)])
    with pytest.raises(ValueError, match="Conflicting"):
        reporting.combine_cells({}, [cell(), cell(identity="different")])


def test_unavailable_original_cell_can_have_explicit_verified_supplement():
    unavailable = {**cell(), "status": "unavailable", "scores": {}}
    result = reporting.combine_cells({}, [unavailable, cell()])["membership:control_direct"]
    assert result["scores"] == {ENDPOINTS[-1]: .8}
    assert result["cells"][0]["status"] == "unavailable"


@pytest.mark.parametrize("kind,flag", [("full_search", "inputs_sources_runtime_reverified"),
                                      ("membership", "inputs_reverified_unchanged"),
                                      ("supplement", "inputs_reverified_unchanged")])
def test_producer_postflight_is_required(kind, flag, tmp_path):
    with pytest.raises(ValueError, match="postflight"):
        reporting.audit_cell(tmp_path, kind, "OrthoBench", {"status": "scored_development_result", flag: False}, {})


def fixture_replay(tmp_path):
    fasta = tmp_path / "input"
    fasta.mkdir()
    for taxon, gene in (("A", "a"), ("B", "b")):
        (fasta / f"{taxon}.faa").write_text(f">{gene}\nACDE\n")
    native = tmp_path / "native"
    phylo, working = native / "orthohmm_phylogeny", native / "orthohmm_working_res"
    phylo.mkdir(parents=True)
    working.mkdir()
    candidate = tmp_path / "candidate.txt"
    candidate.write_text("a b\n")
    (working / "orthohmm_edges_clustered.txt").write_text("a b\n")
    provenance = {"input_proteomes": [{"filename": f"{taxon}.faa", "taxon": taxon,
        "sha256": record(fasta / f"{taxon}.faa")["sha256"]} for taxon in ("A", "B")],
        "root_duplication_rule": "species_overlap", "pair_orthology_rule": "positive_paralogy",
        "species_tree_mode": "infer", "species_tree_rooting": "min_variance"}
    (phylo / "provenance_manifest.json").write_text(json.dumps(provenance))
    (phylo / "reconciliation_summary.json").write_text(json.dumps({"ortholog_pairs": 1, "root_hogs": 1, "candidate_families": 1}))
    (phylo / "species_tree.rooted.nwk").write_text("(A:1,B:1);\n")
    (phylo / reporting.PRODUCTS[0]).write_text("root_hog\tsource_family\tgenes\nH0\tF0\ta,b\n")
    (phylo / reporting.PRODUCTS[1]).write_text("gene_a\tspecies_a\tgene_b\tspecies_b\na\tA\tb\tB\n")
    (phylo / reporting.PRODUCTS[2]).write_text("gene_a\tspecies_a\tgene_b\tspecies_b\tconfidence\na\tA\tb\tB\thigh\n")
    row = {"outputs": [record(phylo / name) for name in reporting.PRODUCTS],
        "coverage": {"candidate_genes": 2, "prediction_genes": 2, "all_candidate_genes_preserved": True}}
    case = {"command": ["python", "replay.py", "--output-directory", str(native), "--fasta-directory", str(fasta),
        "--candidate-clusters", str(candidate), "--root-rule", "species_overlap", "--pair-rule", "positive_paralogy",
        "--species-tree-mode", "infer", "--species-tree-rooting", "min_variance"]}
    return row, case, phylo


def test_native_replay_readback_checks_universe_pair_confidence_and_summary(tmp_path):
    row, case, phylo = fixture_replay(tmp_path)
    result = reporting.replay_native(row, case)
    assert (result["genes"], result["pairs"], result["groups"], result["species"]) == (2, 1, 1, 2)
    (phylo / reporting.PRODUCTS[2]).write_text("gene_a\tspecies_a\tgene_b\tspecies_b\tconfidence\n")
    with pytest.raises(ValueError, match="outputs differ"):
        reporting.replay_native(row, case)


@pytest.mark.parametrize("mutation", ["checksum", "configuration", "tree", "universe", "pair_species", "count", "root"])
def test_replay_corruption_is_not_admitted(tmp_path, mutation):
    row, case, phylo = fixture_replay(tmp_path)
    if mutation in {"checksum", "configuration"}:
        path = phylo / "provenance_manifest.json"
        payload = json.loads(path.read_text())
        if mutation == "checksum":
            payload["input_proteomes"][0]["sha256"] = "foreign"
        else:
            payload["pair_orthology_rule"] = "lca"
        path.write_text(json.dumps(payload))
    elif mutation == "tree":
        (phylo / "species_tree.rooted.nwk").write_text("(A:1,C:1);\n")
    elif mutation == "universe":
        (phylo.parent / "orthohmm_working_res/orthohmm_edges_clustered.txt").write_text("a\n")
    elif mutation == "pair_species":
        (phylo / reporting.PRODUCTS[1]).write_text("gene_a\tspecies_a\tgene_b\tspecies_b\na\tA\tb\tA\n")
    elif mutation == "count":
        (phylo / "reconciliation_summary.json").write_text(json.dumps({"ortholog_pairs": 2, "root_hogs": 1, "candidate_families": 1}))
    else:
        (phylo / reporting.PRODUCTS[0]).write_text("root_hog\tsource_family\tgenes\nH0\tF0\ta\n")
    with pytest.raises(ValueError):
        reporting.replay_native(row, case)


def test_report_keeps_snapshots_failures_and_missing_cells(tmp_path, monkeypatch):
    registry = tmp_path / "registry.json"
    registry.write_text(json.dumps({"rows": [], "limitations": []}))
    source = tmp_path / "membership/OrthoBench"
    source.mkdir(parents=True)
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"source_revision": "fixture"}))
    payload = {"source_revision": "fixture", "execution_plan": record(plan), "arms": [
        {"label": "control_direct", "status": "scored_development_result", "scores": {"OrthoBench F1": .8}},
        {"label": "control_path", "status": "failed", "scores": {"OrthoBench F1": .9}},
        {"label": "broader_direct", "status": "scored_development_result", "scores": {"OrthoBench F1": .85}},
    ]}
    raw_path = source / "results.json"
    raw_path.write_text(json.dumps(payload))
    original = raw_path.read_bytes()
    def audit(repo, kind, dataset, row, execution):
        if row["label"] == "broader_direct":
            raise ValueError("checksum changed")
        return {"status": "independent_content_readback_verified", "scores": {ENDPOINTS[0]: .8}, "products": {}}
    monkeypatch.setattr(reporting, "audit_cell", audit)
    output = tmp_path / "report"
    result = reporting.build_report(tmp_path, registry, source.parent, tmp_path / "missing_supplement.json",
                                   tmp_path / "missing_full", output)
    methods = {m["key"]: m for m in result["methods"]}
    assert methods["membership:control_direct"]["scores"] == {ENDPOINTS[0]: .8}
    assert methods["membership:control_path"]["scores"] == {}
    assert methods["membership:broader_direct"]["scores"] == {}
    assert result["readback_failures"][0]["error"] == "checksum changed"
    assert "NA" in (output / "scores.tsv").read_text()
    assert raw_path.read_bytes() == original
    assert (output / "snapshots/membership_OrthoBench.json").read_bytes() == original
    with pytest.raises(FileExistsError):
        reporting.build_report(tmp_path, registry, source.parent, tmp_path / "missing_supplement.json",
                               tmp_path / "missing_full", output)

import json

import pytest

from benchmark_tools.verify_frozen_phylogeny_install import fixture, validate


def test_fixed_fixture_contains_markers_and_duplicates():
    first = fixture()
    assert first == fixture()
    assert set(first) == {"S1", "S2", "S3", "S4"}
    assert sum(map(len, first.values())) == 16
    for taxon, genes in first.items():
        assert all(len(sequence) == 200 for sequence in genes.values())
        assert genes[f"{taxon}_family2"] == genes[f"{taxon}_family2_duplicate"]
        assert len(set(genes.values())) == 3


def prepared(tmp_path):
    phylo = tmp_path / "phylogeny"
    phylo.mkdir()
    manifest = dict(species_tree_mode="infer", species_tree_rooting="min_variance",
        root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy",
        species_tree_taxa=["S1", "S2"])
    summary = dict(species_tree_families=1, reconciled_families=1, checkpoint_hits=0,
                   species_tree_checkpoint_hit=False, ortholog_pairs=1)
    (phylo / "provenance_manifest.json").write_text(json.dumps(manifest))
    (phylo / "reconciliation_summary.json").write_text(json.dumps(summary))
    (tmp_path / "orthohmm_orthogroups.txt").write_text("OG0: a b\n")
    (phylo / "orthohmm_pairwise_orthologs.tsv").write_text("gene_a\tspecies_a\tgene_b\tspecies_b\na\tS1\tb\tS2\n")
    (phylo / "gene_trees").mkdir()
    (phylo / "gene_trees/CF0.reconciled.nwk").write_text("(a,b);\n")
    return phylo, {"S1": {"a": "AA"}, "S2": {"b": "AA"}}


def test_complete_artifacts(tmp_path):
    _, genes = prepared(tmp_path)
    result = validate(tmp_path, genes)
    assert result["genes"] == 2 and result["pairs"] == 1


@pytest.mark.parametrize("change", ["bypass", "no_markers", "cached", "missing_gene", "duplicate_gene", "same_species", "duplicate_pair", "missing_tree"])
def test_incomplete_or_invalid_phylogeny_rejected(tmp_path, change):
    phylo, genes = prepared(tmp_path)
    if change in {"bypass", "no_markers", "cached"}:
        path = phylo / "reconciliation_summary.json"
        summary = json.loads(path.read_text())
        summary[{"bypass": "reconciled_families", "no_markers": "species_tree_families", "cached": "checkpoint_hits"}[change]] = 1 if change == "cached" else 0
        path.write_text(json.dumps(summary))
    elif change in {"missing_gene", "duplicate_gene"}:
        (tmp_path / "orthohmm_orthogroups.txt").write_text("OG0: a\n" if change == "missing_gene" else "OG0: a b b\n")
    elif change == "same_species":
        path = phylo / "orthohmm_pairwise_orthologs.tsv"
        path.write_text(path.read_text().replace("b\tS2", "b\tS1"))
    elif change == "duplicate_pair":
        path = phylo / "orthohmm_pairwise_orthologs.tsv"
        path.write_text(path.read_text() + "a\tS1\tb\tS2\n")
    else:
        (phylo / "gene_trees/CF0.reconciled.nwk").unlink()
    with pytest.raises(ValueError):
        validate(tmp_path, genes)

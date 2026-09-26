import json

import pytest

from benchmark_tools.audit_phylogeny_structure import audit, pair_readback, tree_leaves
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.fixture
def tree_fixture(tmp_path):
    inputs, output = tmp_path / "input", tmp_path / "phylo"
    inputs.mkdir()
    output.mkdir()
    for name in ("checkpoints", "gene_trees"):
        (output / name).mkdir()
    proteomes = []
    for gene, taxon in (("a", "S1"), ("b", "S2")):
        path = inputs / f"{taxon}.fa"
        path.write_text(f">{gene}\nACDE\n")
        proteomes.append(dict(filename=path.name, taxon=taxon, sha256=record(path)["sha256"]))
    species = output / "species_tree.rooted.nwk"
    species.write_text("[&R] (S1:1,S2:1);\n")
    checkpoint = dict(status="complete", family_id="F", genes=["a", "b"],
                      species_tree_sha256=record(species)["sha256"])
    for suffix, key, tree in (("raw", "raw_tree_sha256", "(g00000000:1,g00000001:1);"),
                             ("rooted", "rooted_tree_sha256", "[&R] (a:1,b:1);"),
                             ("reconciled", "annotated_tree_sha256", "[&R] (a:1,b:1)G0;")):
        path = output / "gene_trees" / f"F.{suffix}.nwk"
        path.write_text(tree + "\n")
        checkpoint[key] = record(path)["sha256"]
    (output / "checkpoints/F.json").write_text(json.dumps(checkpoint))
    summary = dict(root_hogs=1, candidate_families=1, reconciled_families=1,
                   bypassed_families=0, ortholog_pairs=1)
    manifest = dict(results=summary, input_proteomes=proteomes, species_tree_taxa=["S1", "S2"],
                    species_tree_sha256=record(species)["sha256"])
    (output / "provenance_manifest.json").write_text(json.dumps(manifest))
    (output / "reconciliation_summary.json").write_text(json.dumps(summary))
    (output / "orthohmm_root_hogs.tsv").write_text("root_hog\tsource_family\tgenes\nH\tF\ta,b\n")
    (output / "orthohmm_pairwise_orthologs.tsv").write_text(
        "gene_a\tspecies_a\tgene_b\tspecies_b\na\tS1\tb\tS2\n")
    (output / "orthohmm_pairwise_orthologs_confidence.tsv").write_text(
        "gene_a\tspecies_a\tgene_b\tspecies_b\tconfidence\na\tS1\tb\tS2\thigh\n")
    return output, inputs


def test_full_structural_readback(tree_fixture):
    result = audit(*tree_fixture)
    assert (result["genes"], result["parsed_trees"], result["ortholog_pairs"]) == (2, 4, 1)
    assert not result["scientific_scores_admitted"]
    assert not result["reconciliation_semantics_recomputed"]


@pytest.mark.parametrize("mutation", ["foreign_tree", "extra_tree", "missing_tree",
                                    "family", "input", "summary", "species"])
def test_structural_mutations(tree_fixture, mutation):
    output, inputs = tree_fixture
    if mutation == "foreign_tree":
        path = output / "gene_trees/F.rooted.nwk"
        path.write_text("(a:1,foreign:1);\n")
        cp = output / "checkpoints/F.json"
        data = json.loads(cp.read_text())
        data["rooted_tree_sha256"] = record(path)["sha256"]
        cp.write_text(json.dumps(data))
    elif mutation == "extra_tree":
        (output / "gene_trees/extra.nwk").write_text("(a,b);\n")
    elif mutation == "missing_tree":
        (output / "gene_trees/F.raw.nwk").unlink()
    elif mutation == "family":
        path = output / "orthohmm_root_hogs.tsv"
        path.write_text(path.read_text().replace("\tF\t", "\tforeign\t"))
    elif mutation == "input":
        (inputs / "S1.fa").write_text(">a\nAAAA\n")
    elif mutation == "summary":
        path = output / "reconciliation_summary.json"
        data = json.loads(path.read_text())
        data["ortholog_pairs"] = 2
        path.write_text(json.dumps(data))
    else:
        path = output / "provenance_manifest.json"
        data = json.loads(path.read_text())
        data["input_proteomes"][0]["taxon"] = "foreign"
        path.write_text(json.dumps(data))
    with pytest.raises((ValueError, FileNotFoundError)):
        audit(output, inputs)


@pytest.mark.parametrize("body", ["a\tS1\tb\tS2\nhigh", "a\tS1\tb\tS2\n" * 2,
    "b\tS2\ta\tS1\n", "a\twrong\tb\tS2\n", "a\tS1\tx\tS2\n", ""])
def test_bad_pairs(tree_fixture, body):
    output, _ = tree_fixture
    path = output / "orthohmm_pairwise_orthologs.tsv"
    path.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\n" + body)
    with pytest.raises(ValueError):
        pair_readback(output, {"a": "S1", "b": "S2"}, {"a": "F", "b": "F"}, 1)


def test_multiple_trees_rejected(tmp_path):
    path = tmp_path / "tree.nwk"
    path.write_text("(a,b);\n(a,b);\n")
    with pytest.raises(ValueError, match="one tree"):
        tree_leaves(path, {"a", "b"})


@pytest.mark.parametrize("mutation", ["reverse", "duplicate", "foreign", "species",
                                    "same_species", "cross_family", "confidence", "count"])
def test_matching_but_invalid_pair_tables(tree_fixture, mutation):
    output, _ = tree_fixture
    gene_species, gene_family = {"a": "S1", "b": "S2"}, {"a": "F", "b": "F"}
    values, confidence, copies, count = ["a", "S1", "b", "S2"], "high", 1, 1
    if mutation == "reverse":
        values = ["b", "S2", "a", "S1"]
    elif mutation == "duplicate":
        copies = count = 2
    elif mutation == "foreign":
        values[2] = "x"
    elif mutation == "species":
        values[1] = "wrong"
    elif mutation == "same_species":
        gene_species["b"] = values[3] = "S1"
    elif mutation == "cross_family":
        gene_family["b"] = "G"
    elif mutation == "confidence":
        confidence = "invalid"
    else:
        count = 2
    for suffix, extra in (("", []), ("_confidence", [confidence])):
        header = ["gene_a", "species_a", "gene_b", "species_b"] + (["confidence"] if extra else [])
        (output / f"orthohmm_pairwise_orthologs{suffix}.tsv").write_text(
            "\t".join(header) + "\n" + ("\t".join(values + extra) + "\n") * copies)
    with pytest.raises(ValueError):
        pair_readback(output, gene_species, gene_family, count)

from pathlib import Path

import pytest

from orthohmm.phylogeny import (
    PhylogenyConfig,
    PhylogenyConfigurationError,
    canonical_species_name,
    canonical_species_names,
    parse_species_tree,
    parse_gene_tree,
    reconcile_gene_tree,
    root_gene_tree_min_duplication_loss,
    root_tree_min_variance,
    validate_phylogeny_config,
)


pytest.importorskip("dendropy")


PROTEOMES = ["A.faa", "B.fasta", "C.pep"]


def write_tree(tmp_path: Path, newick: str) -> Path:
    path = tmp_path / "species_tree.nwk"
    path.write_text(newick, encoding="utf-8")
    return path


def test_canonical_species_names_strip_supported_extensions():
    assert canonical_species_name("path/Species_name.faa") == "Species_name"
    assert canonical_species_names(PROTEOMES) == ("A", "B", "C")


def test_canonical_species_names_reject_collisions():
    with pytest.raises(PhylogenyConfigurationError, match="unique taxon names"):
        canonical_species_names(["A.fa", "A.faa"])


def test_parse_rooted_tree_and_canonicalize_labels(tmp_path):
    result = parse_species_tree(
        write_tree(tmp_path, "[&R] ((A.faa,B.fasta),C.pep);"),
        PROTEOMES,
    )
    assert result.taxa == ("A", "B", "C")
    assert result.tree.is_rooted


def test_parse_common_unmarked_rooted_tree(tmp_path):
    result = parse_species_tree(write_tree(tmp_path, "((A,B),C);"), PROTEOMES)
    assert result.taxa == ("A", "B", "C")


def test_parse_explicitly_rooted_polytomy(tmp_path):
    result = parse_species_tree(
        write_tree(tmp_path, "[&R] (A,B,C);"),
        PROTEOMES,
    )
    assert result.tree.seed_node.num_child_nodes() == 3


@pytest.mark.parametrize(
    "newick",
    [
        "[&U] ((A,B),C);",
        "(A,B,C);",
    ],
)
def test_reject_unrooted_or_root_ambiguous_tree(tmp_path, newick):
    with pytest.raises(PhylogenyConfigurationError, match="must be rooted"):
        parse_species_tree(write_tree(tmp_path, newick), PROTEOMES)


def test_reject_invalid_newick(tmp_path):
    with pytest.raises(PhylogenyConfigurationError, match="Invalid Newick"):
        parse_species_tree(write_tree(tmp_path, "((A,B),C"), PROTEOMES)


def test_reject_duplicate_tree_taxa(tmp_path):
    with pytest.raises(PhylogenyConfigurationError, match="not unique"):
        parse_species_tree(write_tree(tmp_path, "((A,A),C);"), PROTEOMES)


def test_report_missing_and_extra_taxa(tmp_path):
    with pytest.raises(PhylogenyConfigurationError) as error:
        parse_species_tree(write_tree(tmp_path, "((A,B),D);"), PROTEOMES)
    message = str(error.value)
    assert "missing from tree: C" in message
    assert "not present in input proteomes: D" in message


def test_off_mode_does_not_resolve_dependencies_or_tools():
    config = PhylogenyConfig(
        mode="off",
        aligner="definitely-missing-aligner",
        tree_builder="definitely-missing-tree-builder",
    )
    assert validate_phylogeny_config(config, PROTEOMES) is config


def test_enabled_mode_reports_missing_executable(tmp_path):
    config = PhylogenyConfig(
        mode="reconcile",
        species_tree_mode="supplied",
        species_tree=str(write_tree(tmp_path, "((A,B),C);")),
        aligner="definitely-missing-aligner",
        tree_builder="definitely-missing-tree-builder",
    )
    with pytest.raises(PhylogenyConfigurationError, match="aligner executable"):
        validate_phylogeny_config(config, PROTEOMES)


def test_infer_mode_rejects_supplied_tree(tmp_path):
    config = PhylogenyConfig(
        mode="reconcile",
        species_tree_mode="infer",
        species_tree=str(write_tree(tmp_path, "((A,B),C);")),
    )
    with pytest.raises(PhylogenyConfigurationError, match="cannot be used"):
        validate_phylogeny_config(config, PROTEOMES)


def species_tree(tmp_path):
    return parse_species_tree(
        write_tree(tmp_path, "((A,B),C);"),
        PROTEOMES,
    ).tree


def test_reconcile_labels_speciation_nodes_and_ortholog_pairs(tmp_path):
    gene_tree = parse_gene_tree("((a,b),c);")
    result = reconcile_gene_tree(
        gene_tree,
        species_tree(tmp_path),
        {"a": "A.faa", "b": "B.fasta", "c": "C.pep"},
        family_id="OG0",
    )
    assert result.duplication_count == 0
    assert result.speciation_count == 2
    assert result.root_groups == (("a", "b", "c"),)
    assert result.ortholog_pairs == (("a", "b"), ("a", "c"), ("b", "c"))
    assert result.paralog_pairs == ()
    assert "|S@" in result.annotated_newick


def test_ancient_duplication_splits_root_lineages(tmp_path):
    gene_tree = parse_gene_tree("(((a1,b1),c1),((a2,b2),c2));")
    gene_species = {
        "a1": "A",
        "b1": "B",
        "c1": "C",
        "a2": "A",
        "b2": "B",
        "c2": "C",
    }
    result = reconcile_gene_tree(
        gene_tree,
        species_tree(tmp_path),
        gene_species,
        family_id="OG1",
    )
    assert result.duplication_count == 1
    assert result.root_groups == (
        ("a1", "b1", "c1"),
        ("a2", "b2", "c2"),
    )
    assert ("a1", "b1") in result.ortholog_pairs
    assert ("a1", "b2") in result.paralog_pairs
    assert "|D@S0000" in result.annotated_newick


def test_recent_duplication_does_not_split_root_orthogroup(tmp_path):
    gene_tree = parse_gene_tree("(((a1,a2),b),c);")
    result = reconcile_gene_tree(
        gene_tree,
        species_tree(tmp_path),
        {"a1": "A", "a2": "A", "b": "B", "c": "C"},
    )
    assert result.duplication_count == 1
    assert result.root_groups == (("a1", "a2", "b", "c"),)
    assert ("a1", "a2") not in result.ortholog_pairs
    assert ("a1", "b") in result.ortholog_pairs
    assert ("a2", "b") in result.ortholog_pairs


def test_sparse_root_duplication_is_not_split_without_two_root_lineages(tmp_path):
    gene_tree = parse_gene_tree("(a,((b,c),c2));")
    result = reconcile_gene_tree(
        gene_tree,
        species_tree(tmp_path),
        {"a": "A", "b": "B", "c": "C", "c2": "C"},
    )
    assert result.duplication_count >= 1
    assert result.root_groups == (("a", "b", "c", "c2"),)


def test_differential_loss_preserves_two_ancient_lineages(tmp_path):
    gene_tree = parse_gene_tree("((a1,c1),(b2,c2));")
    result = reconcile_gene_tree(
        gene_tree,
        species_tree(tmp_path),
        {"a1": "A", "c1": "C", "b2": "B", "c2": "C"},
    )
    assert result.root_groups == (("a1", "c1"), ("b2", "c2"))
    assert result.duplication_count == 1


def test_root_duplication_below_speciation_splits_root_hogs(tmp_path):
    gene_tree = parse_gene_tree("(a,((b1,c1),(b2,c2)));")
    result = reconcile_gene_tree(
        gene_tree,
        species_tree(tmp_path),
        {"a": "A", "b1": "B", "c1": "C", "b2": "B", "c2": "C"},
    )
    assert result.root_groups == (
        ("a",),
        ("b1", "c1"),
        ("b2", "c2"),
    )


def test_root_duplication_rule_validation(tmp_path):
    with pytest.raises(ValueError, match="root duplication rule"):
        reconcile_gene_tree(
            parse_gene_tree("(a,b);"),
            species_tree(tmp_path),
            {"a": "A", "b": "B"},
            root_duplication_rule="unknown",
        )


def test_pair_orthology_rule_validation(tmp_path):
    with pytest.raises(ValueError, match="pair orthology rule"):
        reconcile_gene_tree(
            parse_gene_tree("(a,b);"),
            species_tree(tmp_path),
            {"a": "A", "b": "B"},
            pair_orthology_rule="unknown",
        )


def test_positive_paralogy_preserves_pairs_at_mapping_only_conflict(tmp_path):
    mapping = {"a": "A", "b": "B", "c": "C"}
    legacy = reconcile_gene_tree(
        parse_gene_tree("(a,(b,c));"),
        species_tree(tmp_path),
        mapping,
    )
    confidence_aware = reconcile_gene_tree(
        parse_gene_tree("(a,(b,c));"),
        species_tree(tmp_path),
        mapping,
        pair_orthology_rule="positive_paralogy",
    )

    assert ("a", "b") not in legacy.ortholog_pairs
    assert ("a", "b") in confidence_aware.ortholog_pairs
    assert ("a", "b", "medium") in (
        confidence_aware.ortholog_pair_confidence
    )
    assert confidence_aware.uncertain_count == 1
    assert any(
        node.pair_event == "uncertain"
        and node.mapping_conflict
        and node.species_overlap_count == 0
        for node in confidence_aware.nodes
    )


def test_positive_paralogy_still_removes_species_overlap_pairs(tmp_path):
    result = reconcile_gene_tree(
        parse_gene_tree("((a1,b1),(a2,c1));"),
        species_tree(tmp_path),
        {"a1": "A", "b1": "B", "a2": "A", "c1": "C"},
        pair_orthology_rule="positive_paralogy",
    )

    assert ("a1", "c1") in result.paralog_pairs
    assert ("a2", "b1") in result.paralog_pairs
    assert ("a1", "b1") in result.ortholog_pairs
    assert result.uncertain_count == 0


def test_species_overlap_rule_exposes_sparse_duplication_tradeoff(tmp_path):
    result = reconcile_gene_tree(
        parse_gene_tree("(a,((b,c),c2));"),
        species_tree(tmp_path),
        {"a": "A", "b": "B", "c": "C", "c2": "C"},
        root_duplication_rule="species_overlap",
    )

    assert result.root_groups == (("a",), ("b", "c"), ("c2",))


@pytest.mark.parametrize(
    ("support", "expected"),
    [
        ("0.50", (("a", "b", "c", "c2"),)),
        ("95", (("a",), ("b", "c"), ("c2",))),
    ],
)
def test_confidence_root_rule_requires_strong_sparse_split_support(
    tmp_path, support, expected
):
    result = reconcile_gene_tree(
        parse_gene_tree(f"(a,((b,c),c2){support});"),
        species_tree(tmp_path),
        {"a": "A", "b": "B", "c": "C", "c2": "C"},
        root_duplication_rule="confidence",
    )

    assert result.root_groups == expected
    supported_node = next(
        node for node in result.nodes if set(node.genes) == {"b", "c", "c2"}
    )
    assert supported_node.branch_support == pytest.approx(float(support) / (
        100.0 if float(support) > 1.0 else 1.0
    ))


def test_gene_tree_rejects_unknown_gene(tmp_path):
    with pytest.raises(PhylogenyConfigurationError, match="absent"):
        reconcile_gene_tree(
            parse_gene_tree("(a,unknown);"),
            species_tree(tmp_path),
            {"a": "A"},
        )


def test_minimum_duplication_loss_roots_ancient_duplication(tmp_path):
    tree = parse_gene_tree("(((a1,b1),c1),((a2,b2),c2));")
    tree.is_rooted = False
    gene_species = {
        "a1": "A", "b1": "B", "c1": "C",
        "a2": "A", "b2": "B", "c2": "C",
    }
    rooted = root_gene_tree_min_duplication_loss(
        tree, species_tree(tmp_path), gene_species
    )
    result = reconcile_gene_tree(rooted, species_tree(tmp_path), gene_species)
    assert result.duplication_count == 1
    assert result.root_groups == (
        ("a1", "b1", "c1"),
        ("a2", "b2", "c2"),
    )


def test_minimum_duplication_loss_keeps_recent_duplication_in_root_hog(tmp_path):
    tree = parse_gene_tree("(((a1,a2),b),c);")
    tree.is_rooted = False
    gene_species = {"a1": "A", "a2": "A", "b": "B", "c": "C"}
    rooted = root_gene_tree_min_duplication_loss(
        tree, species_tree(tmp_path), gene_species
    )
    result = reconcile_gene_tree(rooted, species_tree(tmp_path), gene_species)
    assert result.duplication_count == 1
    assert result.root_groups == (("a1", "a2", "b", "c"),)


def test_minimum_variance_root_is_deterministic_and_uses_internal_edge():
    dendropy = pytest.importorskip("dendropy")
    tree = dendropy.Tree.get(
        data="((a:1,b:1):3,(c:1,d:1):1,e:1);",
        schema="newick",
        preserve_underscores=True,
        rooting="default-unrooted",
    )

    rooted = root_tree_min_variance(tree)

    root_sides = {
        frozenset(str(leaf.taxon.label) for leaf in child.leaf_iter())
        for child in rooted.seed_node.child_node_iter()
    }
    assert root_sides == {
        frozenset(("a", "b")),
        frozenset(("c", "d", "e")),
    }
    assert rooted.is_rooted

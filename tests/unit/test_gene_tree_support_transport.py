from dataclasses import replace
from pathlib import Path

import pytest

from orthohmm.phylogeny import (
    PhylogenyConfig, PhylogenyConfigurationError, capture_gene_branch_support,
    parse_gene_tree, restore_gene_branch_support, root_gene_tree_min_duplication_loss,
    validate_gene_tree_support_mode, validate_phylogeny_config,
    _gene_tree_split_keys, reconcile_gene_tree,
)
from orthohmm.phylogeny_pipeline import _family_tree_input_hash, _root_external_gene_tree


RAW = "((a:0.1,b:0.2)0.31:0.5,c:0.2,(d:0.1,e:0.1)0.72:0.4);"
SPECIES = "(a:1,(b:1,(c:1,(d:1,e:1):1):1):1);"
OWNERS = {g: g for g in "abcde"}


def distances(tree):
    matrix = tree.phylogenetic_distance_matrix()
    taxa = sorted(tree.taxon_namespace, key=lambda t: t.label)
    return {(a.label, b.label): matrix.distance(a, b) for i, a in enumerate(taxa) for b in taxa[i + 1:]}


def test_production_reproducer_preserves_support_topology_lengths_and_existing_pair_results():
    species = parse_gene_tree(SPECIES)
    old, fixed = parse_gene_tree(RAW), parse_gene_tree(RAW)
    snapshot = capture_gene_branch_support(fixed)
    original_lengths = distances(fixed)
    original_splits = set(_gene_tree_split_keys(fixed)[1].values())
    root_gene_tree_min_duplication_loss(old, species, OWNERS)
    root_gene_tree_min_duplication_loss(fixed, species, OWNERS, gene_tree_support_mode="bipartition")
    assert capture_gene_branch_support(old) != snapshot
    assert capture_gene_branch_support(fixed) == snapshot
    assert set(_gene_tree_split_keys(fixed)[1].values()) == original_splits
    assert distances(fixed) == pytest.approx(original_lengths)
    assert distances(fixed) == pytest.approx(distances(old))
    for tree in (old, fixed):
        result = reconcile_gene_tree(tree, species, OWNERS, root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")
        if tree is old:
            reference = result
        else:
            assert result.ortholog_pairs == reference.ortholog_pairs
            assert result.ortholog_pair_confidence == reference.ortholog_pair_confidence
            assert result.root_groups == reference.root_groups


@pytest.mark.parametrize("rooting", ["midpoint", "min_variance"])
def test_other_unrooted_gene_rooting_routes_preserve_support(rooting):
    expected = capture_gene_branch_support(parse_gene_tree(RAW))
    tree = _root_external_gene_tree(RAW, OWNERS, unrooted_rooting=rooting, gene_tree_support_mode="bipartition")
    assert capture_gene_branch_support(tree) == expected


@pytest.mark.parametrize("target", ["a", "internal"])
def test_supported_edge_root_and_pendant_root_restore_original_splits(target):
    tree = parse_gene_tree(RAW)
    snapshot = capture_gene_branch_support(tree)
    node = next(n for n in tree.preorder_node_iter() if (n.is_leaf() and n.taxon.label == target)
                or (target == "internal" and n.label == "0.31"))
    tree.reroot_at_edge(node.edge, length1=.1, length2=.1, update_bipartitions=False)
    restore_gene_branch_support(tree, snapshot)
    assert capture_gene_branch_support(tree) == snapshot
    if target == "internal":
        assert sorted(child.label for child in tree.seed_node.child_node_iter()) == ["0.31", "0.31"]


def test_explicitly_rooted_inputs_keep_existing_semantics_and_labels():
    raw = "[&R] ((a:1,b:1)0.4:1,(c:1,d:1)0.5:1)0.9;"
    old = _root_external_gene_tree(raw, {g:g for g in "abcd"})
    fixed = _root_external_gene_tree(raw, {g:g for g in "abcd"}, gene_tree_support_mode="bipartition")
    assert fixed.as_string(schema="newick") == old.as_string(schema="newick")


@pytest.mark.parametrize("raw", [
    "((a:1,b:1)31:1,c:1,(d:1,e:1)72:1);",
    "((a:1,b:1)nan:1,c:1,(d:1,e:1)inf:1);",
    "((a:1,b:1)AB:1,c:1,(d:1,e:1)DE:1)ROOT;",
    "((a:1,b:1):1,c:1,(d:1,e:1):1);",
    "((a:1,b:1,c:1)0.2:1,d:1,(e:1,f:1)0.8:1);",
])
def test_percentages_invalid_missing_nonnumeric_and_polytomies(raw):
    expected = capture_gene_branch_support(parse_gene_tree(raw))
    names = {leaf.taxon.label:leaf.taxon.label for leaf in parse_gene_tree(raw).leaf_node_iter()}
    tree = _root_external_gene_tree(raw, names, gene_tree_support_mode="bipartition")
    assert capture_gene_branch_support(tree) == expected
    if "ROOT" in raw:
        assert {n.label for n in tree.preorder_node_iter() if n.label} == {"AB", "DE", "ROOT"}


def test_conflicting_original_split_support_is_not_averaged_or_maximized():
    tree = parse_gene_tree("((a:1,b:1)0.3:1,(c:1,d:1)0.8:1);")
    with pytest.raises(PhylogenyConfigurationError, match="Conflicting"):
        capture_gene_branch_support(tree)
    equivalent = parse_gene_tree("((a:1,b:1)30:1,(c:1,d:1)0.3:1);")
    assert list(capture_gene_branch_support(equivalent)[1].values()) == [.3]


def test_support_restore_rejects_leaf_universe_changes():
    tree = parse_gene_tree(RAW)
    snapshot = capture_gene_branch_support(tree)
    next(tree.leaf_node_iter()).taxon.label = "changed"
    with pytest.raises(PhylogenyConfigurationError, match="universe"):
        restore_gene_branch_support(tree, snapshot)


def test_support_restore_rejects_nonnumeric_identifier_collision_without_overwriting_it():
    tree = parse_gene_tree(RAW)
    snapshot = capture_gene_branch_support(tree)
    supported = next(node for node in tree.preorder_node_iter() if node.label == "0.31")
    supported.label = "named-node"
    with pytest.raises(PhylogenyConfigurationError, match="overwrite a nonnumeric"):
        restore_gene_branch_support(tree, snapshot)
    assert supported.label == "named-node"


def test_zero_support_is_preserved_and_out_of_range_support_is_unavailable():
    tree = parse_gene_tree("((a:1,b:1)0:1,c:1,(d:1,e:1)101:1);")
    snapshot = capture_gene_branch_support(tree)
    assert list(snapshot[1].values()) == [0.]
    fixed = _root_external_gene_tree(tree.as_string(schema="newick"), OWNERS, gene_tree_support_mode="bipartition")
    assert capture_gene_branch_support(fixed) == snapshot


@pytest.mark.parametrize("mode", [None, True, [], "unsupported"])
def test_invalid_modes_fail_before_parsing_external_trees(mode):
    with pytest.raises(PhylogenyConfigurationError, match="mapping mode"):
        _root_external_gene_tree("not newick", {}, gene_tree_support_mode=mode)
    with pytest.raises(PhylogenyConfigurationError):
        validate_gene_tree_support_mode(mode)


def test_opt_in_requires_enabled_phylogeny_and_does_not_invalidate_raw_tree_inputs():
    assert PhylogenyConfig().gene_tree_support_mode == "legacy"
    with pytest.raises(PhylogenyConfigurationError, match="requires phylogeny"):
        validate_phylogeny_config(PhylogenyConfig(gene_tree_support_mode="bipartition"), ())
    config = PhylogenyConfig()
    sequences = {"a":"ACGT", "b":"AGGT"}
    assert _family_tree_input_hash("f", ("a","b"), sequences, config) == _family_tree_input_hash(
        "f", ("a","b"), sequences, replace(config, gene_tree_support_mode="bipartition"))


def test_production_and_replay_parser_modes_are_opt_in():
    from orthohmm.parser import create_parser
    from benchmark_tools.replay_phylogeny import build_parser
    from benchmark_tools.benchmark_production import parse_args
    parser = create_parser()
    assert parser.parse_args(["inputs"]).gene_tree_support_mode == "legacy"
    assert parser.parse_args(["inputs", "--gene_tree_support_mode", "bipartition"]).gene_tree_support_mode == "bipartition"
    replay = build_parser().parse_args(["--fasta-directory", "inputs", "--candidate-clusters", "groups",
                                       "--output-directory", "output", "--json", "result.json",
                                       "--gene-tree-support-mode", "bipartition"])
    assert replay.gene_tree_support_mode == "bipartition"
    args = ["inputs", "output", "result.json", "--cpu", "1"]
    assert parse_args(args).gene_tree_support_mode == "legacy"
    assert parse_args(args + ["--phylogeny", "reconcile", "--gene-tree-support-mode", "bipartition"]).gene_tree_support_mode == "bipartition"
    with pytest.raises(SystemExit):
        parse_args(args + ["--gene-tree-support-mode", "bipartition"])


def test_enabled_validation_preserves_opt_in_support_mode(tmp_path):
    from tests.unit.test_phylogeny_pipeline import _make_fixture
    inputs, _, _, _, config = _make_fixture(tmp_path, "validate_support")
    files = sorted(p.name for p in inputs.iterdir())
    validated = validate_phylogeny_config(replace(config, gene_tree_support_mode="bipartition"), files)
    assert validated.gene_tree_support_mode == "bipartition"


def test_support_mode_is_recorded_and_raw_checkpoint_is_reused(tmp_path):
    from tests.unit.test_phylogeny_pipeline import _make_fixture
    from orthohmm.phylogeny_pipeline import run_phylogeny_stage
    import json
    inputs, output, cluster_path, candidates, config = _make_fixture(tmp_path, "support_cache")
    builder = Path(config.tree_builder)
    builder.write_text("#!/usr/bin/env python3\nimport sys\nif len(sys.argv)==1: print('fake-fasttree 1.0')\n"
        "else: print('((g00000000:1,g00000002:1)0.2:1,(g00000001:1,g00000003:1)0.8:1,(g00000004:1,g00000005:1)0.6:1);')\n")
    config = replace(config, root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")
    files = sorted(p.name for p in inputs.iterdir())
    run_phylogeny_stage(str(inputs), str(output), files, config, 1)
    phylo = output / "orthohmm_phylogeny"
    products = {name:(phylo/name).read_bytes() for name in ("orthohmm_root_hogs.tsv", "orthohmm_pairwise_orthologs.tsv", "orthohmm_pairwise_orthologs_confidence.tsv")}
    raw = (phylo/"gene_trees/Family0000000.raw.nwk").read_bytes()
    cluster_path.write_text(candidates)
    fixed = run_phylogeny_stage(str(inputs), str(output), files,
        replace(config, gene_tree_support_mode="bipartition"), 1)
    assert fixed.checkpoint_hits >= 1
    assert (phylo/"gene_trees/Family0000000.raw.nwk").read_bytes() == raw
    assert {name:(phylo/name).read_bytes() for name in products} == products
    assert json.loads((phylo/"provenance_manifest.json").read_text())["gene_tree_support_mode"] == "bipartition"
    assert json.loads((phylo/"checkpoints/Family0000000.json").read_text())["gene_tree_support_mode"] == "bipartition"

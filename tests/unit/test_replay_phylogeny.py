from pathlib import Path
import json

import pytest

from benchmark_tools import replay_phylogeny


def test_parser_accepts_checkpoint_source():
    args = replay_phylogeny.build_parser().parse_args([
        "--fasta-directory", "proteomes",
        "--candidate-clusters", "clusters.txt",
        "--output-directory", "output",
        "--json", "result.json",
        "--checkpoint-source", "previous",
        "--root-rule", "species_overlap",
        "--pair-rule", "positive_paralogy",
        "--species-tree-rooting", "min_variance",
        "--membership-constraints", "merges.json",
    ])

    assert args.checkpoint_source == Path("previous")
    assert args.root_rule == "species_overlap"
    assert args.pair_rule == "positive_paralogy"
    assert args.species_tree_rooting == "min_variance"
    assert args.membership_constraints == Path("merges.json")


def test_load_membership_constraints(tmp_path):
    path = tmp_path / "merges.json"
    payload = [{"source_genes": ["a"], "target_genes": ["b"]}]
    path.write_text(json.dumps(payload))

    assert replay_phylogeny.load_membership_constraints(path) == payload


def test_load_membership_constraints_rejects_missing_groups(tmp_path):
    path = tmp_path / "merges.json"
    path.write_text('[{"source_genes": ["a"]}]')

    with pytest.raises(SystemExit, match="target_genes"):
        replay_phylogeny.load_membership_constraints(path)


def test_seed_checkpoint_output_copies_required_tree(tmp_path):
    source = tmp_path / "source" / "orthohmm_phylogeny"
    (source / "checkpoints").mkdir(parents=True)
    (source / "gene_trees").mkdir()
    (source / "checkpoints" / "Family0.json").write_text("{}\n")
    (source / "gene_trees" / "Family0.raw.nwk").write_text("(a,b);\n")
    output = tmp_path / "output"

    selected = replay_phylogeny.seed_checkpoint_output(source.parent, output)

    assert selected == source.resolve()
    assert (
        output / "orthohmm_phylogeny" / "checkpoints" / "Family0.json"
    ).read_text() == "{}\n"
    assert (
        output / "orthohmm_phylogeny" / "gene_trees" / "Family0.raw.nwk"
    ).read_text() == "(a,b);\n"


def test_seed_checkpoint_output_rejects_incomplete_source(tmp_path):
    source = tmp_path / "source"
    source.mkdir()

    with pytest.raises(SystemExit, match="missing required directories"):
        replay_phylogeny.seed_checkpoint_output(source, tmp_path / "output")

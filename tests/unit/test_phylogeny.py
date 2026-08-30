from pathlib import Path

import pytest

from orthohmm.phylogeny import (
    PhylogenyConfig,
    PhylogenyConfigurationError,
    canonical_species_name,
    canonical_species_names,
    parse_species_tree,
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

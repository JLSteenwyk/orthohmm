from pathlib import Path
import json

import pytest

from orthohmm.phylogeny import PhylogenyConfig, PhylogenyConfigurationError
from orthohmm.phylogeny_pipeline import (
    _aligner_command,
    family_requires_reconciliation,
    run_phylogeny_stage,
    select_species_tree_families,
)


def _write_executable(path: Path, source: str) -> str:
    path.write_text(source, encoding="utf-8")
    path.chmod(0o755)
    return str(path)


def _make_fixture(tmp_path: Path, name: str):
    root = tmp_path / name
    inputs = root / "inputs"
    output = root / "output"
    working = output / "orthohmm_working_res"
    inputs.mkdir(parents=True)
    working.mkdir(parents=True)
    proteins = {
        "A.faa": {"a1": "MKTUA", "a2": "MKTAT", "single_a": "MCCCC"},
        "B.faa": {"b1": "MKTGA", "b2": "MKTGT", "single_b": "MCCCD"},
        "C.faa": {"c1": "MKAGA", "c2": "MKAGT", "single_c": "MCCCE"},
    }
    for filename, records in proteins.items():
        lines = [f">{gene}\n{sequence}" for gene, sequence in records.items()]
        (inputs / filename).write_text("\n".join(lines) + "\n", encoding="utf-8")
    candidate_clusters = (
        "a1 a2 b1 b2 c1 c2\n"
        "single_a single_b single_c\n"
    )
    cluster_path = working / "orthohmm_edges_clustered.txt"
    cluster_path.write_text(candidate_clusters, encoding="utf-8")
    species_tree = root / "species.nwk"
    species_tree.write_text("((A,B),C);\n", encoding="utf-8")

    aligner = _write_executable(
        root / "fake_mafft",
        "#!/usr/bin/env python3\n"
        "import pathlib, sys\n"
        "if '--version' in sys.argv:\n"
        "    print('fake-mafft 1.0')\n"
        "else:\n"
        "    print(pathlib.Path(sys.argv[-1]).read_text(), end='')\n",
    )
    tree_builder = _write_executable(
        root / "fake_FastTree",
        "#!/usr/bin/env python3\n"
        "import pathlib, sys\n"
        "if len(sys.argv) == 1:\n"
        "    print('fake-fasttree 1.0')\n"
        "elif pathlib.Path(sys.argv[-1]).name == 'species_tree_alignment.faa':\n"
        "    print('[&R] ((s00000000,s00000001),s00000002);')\n"
        "else:\n"
        "    print('[&R] (((g00000000,g00000002),g00000004),"
        "((g00000001,g00000003),g00000005));')\n",
    )
    config = PhylogenyConfig(
        mode="reconcile",
        species_tree_mode="supplied",
        species_tree=str(species_tree),
        aligner=aligner,
        tree_builder=tree_builder,
    )
    return inputs, output, cluster_path, candidate_clusters, config


def test_family_reconciliation_bypass_rules():
    species = {"a": "A", "b": "B", "c": "C", "a2": "A"}
    assert not family_requires_reconciliation(("a",), species)
    assert not family_requires_reconciliation(("a", "b", "c"), species)
    assert family_requires_reconciliation(("a", "a2", "b"), species)


def test_mafft_backend_accepts_nonstandard_protein_symbols(tmp_path):
    command = _aligner_command("/tools/mafft", tmp_path / "family.faa")
    assert "--anysymbol" in command


def test_species_tree_family_selection_is_ranked_and_requires_all_taxa():
    records = [
        ("Family2", ("a2", "b2")),
        ("Family1", ("a1", "b1", "c1")),
    ]
    species = {
        "a1": "A", "b1": "B", "c1": "C", "a2": "A", "b2": "B"
    }
    sequences = {gene: "M" * 10 for gene in species}
    selected = select_species_tree_families(
        records, species, sequences, ("A", "B", "C"), max_families=1
    )
    assert selected == [("Family1", ("a1", "b1", "c1"))]
    with pytest.raises(PhylogenyConfigurationError, match="represent these taxa"):
        select_species_tree_families(
            [("Family2", ("a2", "b2"))],
            species,
            sequences,
            ("A", "B", "C"),
        )


def test_species_tree_selection_adds_sparse_taxon_placement_markers():
    records = [
        ("Core", ("a1", "b1", "c1")),
        ("PlaceD", ("a2", "d1")),
        ("PlaceE", ("b2", "e1")),
        ("PlaceF", ("c2", "f1")),
    ]
    species = {
        "a1": "A", "a2": "A", "b1": "B", "b2": "B",
        "c1": "C", "c2": "C", "d1": "D", "e1": "E", "f1": "F",
    }
    sequences = {gene: "M" * 10 for gene in species}

    selected = select_species_tree_families(
        records,
        species,
        sequences,
        ("A", "B", "C", "D", "E", "F"),
        max_families=1,
    )

    assert selected == records


def test_species_tree_selection_rejects_disconnected_sparse_markers():
    records = [
        ("Core", ("a1", "b1", "c1")),
        ("Detached", ("d1", "e1")),
    ]
    species = {gene: gene[0].upper() for _, genes in records for gene in genes}
    sequences = {gene: "M" * 10 for gene in species}

    with pytest.raises(PhylogenyConfigurationError, match="no connected"):
        select_species_tree_families(
            records,
            species,
            sequences,
            ("A", "B", "C", "D", "E"),
        )


def test_supplied_tree_stage_end_to_end_and_checkpoint_restart(tmp_path):
    inputs, output, cluster_path, candidates, config = _make_fixture(
        tmp_path, "first"
    )
    files = ["A.faa", "B.faa", "C.faa"]
    first = run_phylogeny_stage(
        str(inputs), str(output), files, config, cpu=2
    )
    assert first.candidate_families == 2
    assert first.reconciled_families == 1
    assert first.bypassed_families == 1
    assert first.root_hogs == 3
    assert first.duplications == 1
    assert first.checkpoint_hits == 0
    assert first.remapped_checkpoint_hits == 0
    assert cluster_path.read_text(encoding="utf-8") == (
        "a1 b1 c1\na2 b2 c2\nsingle_a single_b single_c\n"
    )

    pairwise = (
        output / "orthohmm_phylogeny" / "orthohmm_pairwise_orthologs.tsv"
    ).read_text(encoding="utf-8")
    assert "a1\tA\tb1\tB" in pairwise
    assert "a1\tA\tb2\tB" not in pairwise
    assert "single_a\tA\tsingle_b\tB" in pairwise
    confidence_pairwise = (
        output
        / "orthohmm_phylogeny"
        / "orthohmm_pairwise_orthologs_confidence.tsv"
    ).read_text(encoding="utf-8")
    assert "gene_a\tspecies_a\tgene_b\tspecies_b\tconfidence" in (
        confidence_pairwise
    )
    assert "a1\tA\tb1\tB\thigh" in confidence_pairwise
    assert (
        output
        / "orthohmm_phylogeny"
        / "gene_trees"
        / "Family0000000.reconciled.nwk"
    ).is_file()

    cluster_path.write_text(candidates, encoding="utf-8")
    second = run_phylogeny_stage(
        str(inputs), str(output), files, config, cpu=2
    )
    assert second.checkpoint_hits == 1
    assert second == first.__class__(**{**first.__dict__, "checkpoint_hits": 1})


def test_gene_tree_checkpoint_survives_candidate_renumbering(tmp_path):
    inputs, output, cluster_path, candidates, config = _make_fixture(
        tmp_path, "renumbered"
    )
    files = ["A.faa", "B.faa", "C.faa"]
    first = run_phylogeny_stage(
        str(inputs), str(output), files, config, cpu=2
    )
    assert first.checkpoint_hits == 0

    Path(config.aligner).write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if '--version' in sys.argv:\n"
        "    print('fake-mafft 1.0')\n"
        "else:\n"
        "    raise SystemExit(91)\n",
        encoding="utf-8",
    )
    Path(config.tree_builder).write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if len(sys.argv) == 1:\n"
        "    print('fake-fasttree 1.0')\n"
        "else:\n"
        "    raise SystemExit(92)\n",
        encoding="utf-8",
    )
    cluster_path.write_text(
        "single_a single_b single_c\n"
        "a1 a2 b1 b2 c1 c2\n",
        encoding="utf-8",
    )

    renumbered = run_phylogeny_stage(
        str(inputs), str(output), files, config, cpu=2
    )

    assert renumbered.checkpoint_hits == 1
    assert renumbered.remapped_checkpoint_hits == 1


def test_gene_tree_checkpoint_survives_species_tree_change(tmp_path):
    inputs, output, cluster_path, candidates, config = _make_fixture(
        tmp_path, "tree-change"
    )
    files = ["A.faa", "B.faa", "C.faa"]
    run_phylogeny_stage(str(inputs), str(output), files, config, cpu=2)

    # Keep version probes working but make repeated alignment/tree inference fail.
    Path(config.aligner).write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if '--version' in sys.argv:\n"
        "    print('fake-mafft 1.0')\n"
        "else:\n"
        "    raise SystemExit(91)\n",
        encoding="utf-8",
    )
    Path(config.tree_builder).write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if len(sys.argv) == 1:\n"
        "    print('fake-fasttree 1.0')\n"
        "else:\n"
        "    raise SystemExit(92)\n",
        encoding="utf-8",
    )
    Path(config.species_tree).write_text("((A,C),B);\n", encoding="utf-8")
    cluster_path.write_text(candidates, encoding="utf-8")

    rerooted = run_phylogeny_stage(
        str(inputs), str(output), files, config, cpu=2
    )
    assert rerooted.checkpoint_hits == 1


def test_parallel_results_are_deterministic(tmp_path):
    fixture_a = _make_fixture(tmp_path, "run_a")
    fixture_b = _make_fixture(tmp_path, "run_b")
    files = ["A.faa", "B.faa", "C.faa"]
    result_a = run_phylogeny_stage(
        str(fixture_a[0]), str(fixture_a[1]), files, fixture_a[4], cpu=1
    )
    result_b = run_phylogeny_stage(
        str(fixture_b[0]), str(fixture_b[1]), files, fixture_b[4], cpu=4
    )
    assert result_a.__dict__ | {"output_directory": "ignored"} == (
        result_b.__dict__ | {"output_directory": "ignored"}
    )
    for filename in (
        "orthohmm_root_hogs.tsv",
        "orthohmm_pairwise_orthologs.tsv",
        "orthohmm_pairwise_orthologs_confidence.tsv",
        "orthohmm_reconciliation_nodes.tsv",
        "orthohmm_hierarchical_orthogroups.tsv",
    ):
        assert (
            fixture_a[1] / "orthohmm_phylogeny" / filename
        ).read_text(encoding="utf-8") == (
            fixture_b[1] / "orthohmm_phylogeny" / filename
        ).read_text(encoding="utf-8")
    summaries = []
    for fixture in (fixture_a, fixture_b):
        summary = json.loads(
            (fixture[1] / "orthohmm_phylogeny" / "reconciliation_summary.json")
            .read_text(encoding="utf-8")
        )
        summary.pop("output_directory")
        summaries.append(summary)
    assert summaries[0] == summaries[1]


def test_inferred_species_tree_end_to_end_and_restart(tmp_path):
    inputs, output, cluster_path, candidates, supplied_config = _make_fixture(
        tmp_path, "inferred"
    )
    for filename, gene, sequence in (
        ("A.faa", "marker_a", "MNNNN"),
        ("B.faa", "marker_b", "MNNND"),
        ("C.faa", "marker_c", "MNNNE"),
    ):
        with (inputs / filename).open("a", encoding="utf-8") as handle:
            handle.write(f">{gene}\n{sequence}\n")
    candidates += "marker_a marker_b marker_c\n"
    cluster_path.write_text(candidates, encoding="utf-8")
    config = PhylogenyConfig(
        mode="reconcile",
        species_tree_mode="infer",
        aligner=supplied_config.aligner,
        tree_builder=supplied_config.tree_builder,
    )
    files = ["A.faa", "B.faa", "C.faa"]
    first = run_phylogeny_stage(
        str(inputs), str(output), files, config, cpu=3
    )
    assert first.species_tree_families == 2
    assert first.species_tree_checkpoint_hit is False
    assert first.reconciled_families == 1
    assert first.root_hogs == 4
    inferred_tree = output / "orthohmm_phylogeny" / "species_tree.rooted.nwk"
    assert inferred_tree.is_file()
    assert "A" in inferred_tree.read_text(encoding="utf-8")

    cluster_path.write_text(candidates, encoding="utf-8")
    second = run_phylogeny_stage(
        str(inputs), str(output), files, config, cpu=1
    )
    assert second.species_tree_checkpoint_hit is True
    assert second.checkpoint_hits == 1

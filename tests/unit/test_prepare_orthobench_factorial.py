import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools.prepare_orthobench_factorial import indexed_species, plan_cells, prepare_partition
from benchmark_tools.replay_phylogeny import load_membership_constraints


def test_eight_cells_have_own_constraints_and_no_reference_arguments():
    cells = plan_cells(Path("/out"), Path("/input"), Path("/launcher"), 32)
    assert len(cells) == 8 and len({c["label"] for c in cells}) == 8
    assert len({(c["profile_expansion"], c["candidate_expansion"], c["reconciliation"]) for c in cells}) == 8
    for row in cells:
        if not row["reconciliation"]:
            assert "argv" not in row and row["prediction"] == row["candidate_partition"]
            continue
        argv = row["argv"]
        assert "--official-benchmark" not in argv
        assert argv[argv.index("--species-tree-mode") + 1] == "infer"
        assert argv[argv.index("--species-tree-rooting") + 1] == "min_variance"
        assert argv[argv.index("--root-rule") + 1] == "species_overlap"
        assert ("--membership-constraints" in argv) == row["candidate_expansion"]
        if row["candidate_expansion"]:
            assert f"p{int(row['profile_expansion'])}_c1" in argv[argv.index("--membership-constraints") + 1]


def test_cache_species_labels_must_match_fasta_classes_bijectively():
    payload = {"all_gene_ids": ["c", "a", "b"], "gene_to_species": {"a": "0", "b": "0", "c": "1"}}
    owners = {"a": "a.fa", "b": "a.fa", "c": "b.fa"}
    names, species = indexed_species(payload, owners)
    assert names == ["a", "b", "c"] and species.tolist() == [0, 0, 1]
    payload["gene_to_species"]["c"] = "0"
    with pytest.raises(ValueError, match="species classes"):
        indexed_species(payload, owners)
    payload["all_gene_ids"].append("a")
    with pytest.raises(ValueError, match="gene universe"):
        indexed_species(payload, owners)


def test_each_expansion_rebuilds_constraints_and_preserves_seed(tmp_path):
    names = ["a", "b", "c"]
    seed = tmp_path / "seed.txt"
    seed.write_text("a\nb\nc\n")
    calls = []
    def engine(target, observed_names, species, hits, profile):
        calls.append((target, profile))
        work = Path(target) / "orthohmm_working_res"
        (work / "orthohmm_edges_clustered.txt").write_text("a b\nc\n")
        (work / "phylogeny_candidate_merges.json").write_text(json.dumps([{"source_genes": ["a"], "target_genes": ["b"]}]))
        return {"profile": profile, "_membership_constraints": ["not duplicated in report"]}
    fixed = prepare_partition(seed, tmp_path / "no_expansion", names, np.array([0, 1, 2]), (), False, engine, load_membership_constraints)
    assert calls == [] and "membership_constraints" not in fixed
    for label in ("profile_off", "profile_on"):
        expanded = prepare_partition(seed, tmp_path / label, names, np.array([0, 1, 2]), (), True, engine, load_membership_constraints)
        assert "_membership_constraints" not in expanded["expansion"]
        assert label in expanded["membership_constraints"]["path"]
    assert len(calls) == 2 and all(c[1] == "satellite_v2" for c in calls)
    assert seed.read_text() == "a\nb\nc\n"
    with pytest.raises(FileExistsError):
        prepare_partition(seed, tmp_path / "profile_on", names, np.array([0, 1, 2]), (), True, engine, load_membership_constraints)


def test_incomplete_seed_partition_rejected_before_output_creation(tmp_path):
    seed = tmp_path / "seed.txt"
    seed.write_text("a\n")
    with pytest.raises(ValueError, match="does not cover"):
        prepare_partition(seed, tmp_path / "out", ["a", "b"], np.array([0, 1]), (), False, None, None)
    assert not (tmp_path / "out").exists()

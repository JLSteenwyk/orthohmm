from copy import deepcopy
import csv
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import local_rescue_readback as readback
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from orthohmm.local_rescue import evaluate_local_context
from orthohmm.phylogeny import PhylogenyConfig
from tests.unit.test_local_rescue_profiles import workflow_arguments


@pytest.fixture(scope="module")
def accepted(tmp_path_factory):
    if not shutil.which("mafft") or not shutil.which("FastTree"):
        pytest.skip("Requires real MAFFT/FastTree")
    from benchmark_tools.verify_frozen_phylogeny_install import fixture
    from orthohmm.phylogeny_pipeline import _load_sequence_data, run_phylogeny_stage
    from orthohmm.search.profile_expansion import load_global_sequence_database
    root = tmp_path_factory.mktemp("local_rescue_original")
    fasta = root / "input"
    fasta.mkdir()
    for taxon, genes in fixture().items():
        (fasta / f"{taxon}.faa").write_text("".join(f">{gene}\n{seq}\n" for gene, seq in genes.items()))
    files = sorted(p.name for p in fasta.iterdir())
    sequences, species = _load_sequence_data(str(fasta), files)
    names = sorted(sequences)
    groups = [[g for g in names if g.endswith(f"family{i}")] for i in range(3)]
    groups.append([g for g in names if g.endswith("duplicate")])
    baseline = root / "baseline"
    partition = baseline / "orthohmm_working_res/orthohmm_edges_clustered.txt"
    partition.parent.mkdir(parents=True)
    partition.write_text("".join(" ".join(group) + "\n" for group in groups))
    config = PhylogenyConfig(mode="reconcile", species_tree_mode="infer", species_tree_rooting="min_variance",
        root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy",
        aligner=shutil.which("mafft"), tree_builder=shutil.which("FastTree"))
    run_phylogeny_stage(str(fasta), str(baseline), files, config, 1)
    phylo = baseline / "orthohmm_phylogeny"
    genes = set(groups[2] + groups[3])
    with (phylo / "orthohmm_pairwise_orthologs_confidence.tsv").open() as stream:
        pairs = [(r["gene_a"], r["gene_b"]) for r in csv.DictReader(stream, delimiter="\t")
                 if r["confidence"] == "high" and r["gene_a"] in genes and r["gene_b"] in genes]
    context = {"context_id": "LocalContext000001", "family_a": 2, "family_b": 3,
               "genes_a": groups[2], "genes_b": groups[3]}
    database = load_global_sequence_database(str(fasta), files, names)
    original = root / "original_joint"
    proposal = evaluate_local_context(context, names, database, species, sequences,
        phylo / "species_tree.rooted.nwk", files, config, original, [], pairs)
    assert proposal["accepted"]
    return {"root": root, "baseline": phylo, "proposal": original / "proposal.json",
            "arguments": [context, names, database, species, sequences, record(phylo / "species_tree.rooted.nwk"),
                          files, config, [], pairs], "proposal_object": proposal}


def test_actual_profiles_and_cached_joint_tree_reproduce_without_tool_execution(accepted, tmp_path, monkeypatch):
    def forbidden(*a, **kw):
        pytest.fail("Readback must not infer replacement alignments/trees")
    monkeypatch.setattr("orthohmm.phylogeny_pipeline._run_to_file", forbidden)
    result = readback.read_context(accepted["proposal"], *accepted["arguments"], tmp_path / "readback", cpu=1)
    assert result["accepted"] and result["tree_inferred"] and result["profile_accepted"]
    assert result["joint_admission"] == accepted["proposal_object"]["joint_admission"]
    assert len(result["checked"]) == 8
    assert json.loads(accepted["proposal"].read_text()) == accepted["proposal_object"]
    with pytest.raises(FileExistsError):
        readback.read_context(accepted["proposal"], *accepted["arguments"], tmp_path / "readback", cpu=1)


@pytest.mark.parametrize("kind", ["context", "tree_identity", "constraints", "profile", "decision", "accepted", "tokens", "commands", "unfinished"])
def test_proposal_fields_cannot_override_recomputed_admission(accepted, tmp_path, kind):
    source = tmp_path / "source"
    shutil.copytree(accepted["proposal"].parent, source)
    proposal = deepcopy(accepted["proposal_object"])
    # Commands are source-bound, so explicitly rebind only these copy paths.
    for command in proposal["commands"].values():
        command[-1] = command[-1].replace(str(accepted["proposal"].parent), str(source))
    if kind == "context": proposal["context"]["family_a"] += 1
    elif kind == "tree_identity": proposal["species_tree"]["sha256"] = "wrong"
    elif kind == "constraints": proposal["constraints"] = [{"unexpected": True}]
    elif kind == "profile": proposal["profile_validation"]["hits"][0]["raw_score"] += 1
    elif kind == "decision": proposal["joint_admission"]["novel_cross_high_pairs"] = []
    elif kind == "accepted": proposal["accepted"] = False
    elif kind == "tokens": proposal["token_to_gene"]["g00000000"] = "other"
    elif kind == "commands": proposal["commands"]["gene_tree"].insert(1, "-wag")
    else: proposal["status"] = "joint_tree_running"
    (source / "proposal.json").write_text(json.dumps(proposal))
    with pytest.raises(ValueError):
        readback.read_context(source / "proposal.json", *accepted["arguments"], tmp_path / "readback", cpu=1)


@pytest.mark.parametrize("kind", ["raw", "rooted", "reconciled", "candidate", "alignment", "input_hash", "genes"])
def test_tampered_original_tree_or_sequence_artifacts_fail_before_reconciliation(accepted, tmp_path, monkeypatch, kind):
    source = tmp_path / "source"
    shutil.copytree(accepted["proposal"].parent, source)
    family = "LocalContext000001"
    if kind in ("raw", "rooted", "reconciled"):
        path = source / "gene_trees" / (family + "." + kind + ".nwk")
        path.write_text(path.read_text() + "\n")
    elif kind in ("candidate", "alignment"):
        path = source / ("candidate_fastas" if kind == "candidate" else "alignments") / (family + ".faa")
        path.write_text(path.read_text().replace("g00000000", "unknown_token"))
    else:
        path = source / "checkpoints" / (family + ".json")
        value = json.loads(path.read_text())
        value["tree_input_sha256" if kind == "input_hash" else "genes"] = "wrong"
        path.write_text(json.dumps(value))
    def forbidden(*a, **kw):
        pytest.fail("Unauthenticated source tree must never reach reconciliation")
    monkeypatch.setattr(readback, "_reconcile_family", forbidden)
    with pytest.raises(ValueError):
        readback.read_context(source / "proposal.json", *accepted["arguments"], tmp_path / "readback", cpu=1)


def test_rejected_profiles_are_read_back_without_tree_or_fabricated_acceptance(tmp_path, monkeypatch):
    args = workflow_arguments(tmp_path, unrelated=True)
    args[9] = [{"source_genes": ["a1"], "target_genes": ["a2"], "margin": float("inf")}]
    proposal = evaluate_local_context(*args)
    def forbidden(*a, **kw):
        pytest.fail("Profile rejection must not read or infer a tree")
    monkeypatch.setattr(readback, "_reconcile_family", forbidden)
    result = readback.read_context(args[8] / "proposal.json", *args[:5], record(args[5]),
        args[6], args[7], args[9], args[10], tmp_path / "readback", cpu=1)
    assert result["status"] == "profile_rejection_read_back" and not result["accepted"]
    assert not result["tree_inferred"] and not proposal["tree_inferred"]
    assert args[9][0]["margin"] == float("inf")


def test_composite_cache_contains_only_read_back_admissions_and_preserves_sources(accepted, tmp_path):
    from orthohmm.phylogeny_pipeline import _load_reusable_raw_trees
    original = [record(p) for p in sorted(accepted["baseline"].rglob("*")) if p.is_file()]
    proposal = readback.read_context(accepted["proposal"], *accepted["arguments"], tmp_path / "readback", cpu=1)
    destination = tmp_path / "composite"
    result = readback.build_composite_cache(accepted["baseline"], [proposal], destination)
    assert result["accepted_contexts"] == ["LocalContext000001"]
    assert len(result["files"]) == len(original) + 2
    assert not (destination / "alignments/LocalContext000001.faa").exists()
    assert not (destination / "species_tree_inference/LocalContext000001.faa").exists()
    genes = tuple(proposal["genes"])
    assert _load_reusable_raw_trees(destination)[genes].family_id == "LocalContext000001"
    assert [record(p["path"]) for p in original] == original
    with pytest.raises(FileExistsError):
        readback.build_composite_cache(accepted["baseline"], [proposal], destination)


@pytest.mark.parametrize("kind", ["unchecked", "changed_pin", "changed_genes", "id_collision"])
def test_composite_rejects_unchecked_tampered_or_ambiguous_admissions(accepted, tmp_path, kind):
    proposal = readback.read_context(accepted["proposal"], *accepted["arguments"], tmp_path / "readback", cpu=1)
    baseline = accepted["baseline"]
    if kind == "unchecked": proposal["status"] = "running"
    elif kind == "changed_pin": proposal["checked"][0]["sha256"] = "wrong"
    elif kind == "changed_genes": proposal["genes"].append("outside")
    else:
        baseline = tmp_path / "collision"
        shutil.copytree(accepted["baseline"], baseline)
        (baseline / "gene_trees").mkdir(exist_ok=True)
        (baseline / "gene_trees/LocalContext000001.raw.nwk").write_text("collision")
    with pytest.raises(ValueError):
        readback.build_composite_cache(baseline, [proposal], tmp_path / "composite")


def test_actual_scientific_compatibility_with_original_probe_checkout():
    source = Path(readback.__file__).resolve().parents[1]
    original = source / "benchmarks/work/local_joint_rescue_v2_20261018"
    if not original.is_dir():
        pytest.skip("Requires retained immutable probe checkout")
    actual = readback.scientific_compatibility(original, source)
    assert actual["status"] == "original_local_scientific_code_unchanged" and actual["checked"]


@pytest.mark.parametrize("kind", ["core", "profile_helper", "unexpected_module"])
def test_scientific_changes_cannot_hide_inside_readback_namespace(tmp_path, kind):
    old, new = tmp_path / "old", tmp_path / "new"
    for root in (old, new):
        (root / "orthohmm/search").mkdir(parents=True)
        (root / "benchmark_tools").mkdir()
        (root / "orthohmm/core.py").write_text("value = 1\n")
        (root / "orthohmm/search/profile_expansion.py").write_text("def build(): return 1\ndef expand_profiles(): return 1\n")
        for name in ("run_local_rescue_probe.py", "replay_phylogeny.py"):
            (root / "benchmark_tools" / name).write_text("value = 1\n")
    (new / "orthohmm/search/profile_expansion.py").write_text("def build(): return 1\ndef expand_profiles(mode='off'): return 2\n")
    assert readback.scientific_compatibility(old, new)["status"] == "original_local_scientific_code_unchanged"
    if kind == "core": (new / "orthohmm/core.py").write_text("value = 2\n")
    elif kind == "profile_helper": (new / "orthohmm/search/profile_expansion.py").write_text("def build(): return 2\ndef expand_profiles(): return 1\n")
    else: (new / "orthohmm/search/unknown.py").write_text("value = 1\n")
    with pytest.raises(ValueError, match="scientific code differs"):
        readback.scientific_compatibility(old, new)

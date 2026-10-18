from copy import deepcopy
import csv
from dataclasses import replace
import json
import math
import shutil

import numpy as np
import pytest

from benchmark_tools.hmm_generalization_holdout import species_sequences
from orthohmm.local_rescue import evaluate_local_context, validate_local_profiles, write_effective_candidates
from orthohmm.phylogeny import PhylogenyConfig


def sequence_fixture(unrelated=False):
    rng = np.random.default_rng(18102026)
    alphabet = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    ancestor = rng.choice(alphabet, 200)
    sequences, species = {}, {}
    for index in range(3):
        sequence = ancestor.copy()
        sequence[rng.choice(200, 10, replace=False)] = rng.choice(alphabet, 10)
        for side in "ab":
            gene = f"{side}{index + 1}"
            value = rng.choice(alphabet, 200) if unrelated and side == "b" else sequence.copy()
            if side == "b" and not unrelated:
                value[199] = alphabet[(int(np.where(alphabet == value[199])[0][0]) + 1) % 20]
            sequences[gene], species[gene] = "".join(value), f"sp{index + 1}"
    # Included in significance but never targeted or added to a source profile.
    sequences["outside"] = "W" * 1000
    species["outside"] = "sp1"
    names = sorted(sequences)
    database = species_sequences(names, [sequences[gene] for gene in names], "complete")
    context = {"context_id": "LocalContext000001", "genes_a": ["a1", "a2", "a3"],
               "genes_b": ["b1", "b2", "b3"], "family_a": 0, "family_b": 1}
    return context, names, database, species, sequences


def test_real_reciprocal_profiles_use_full_database_and_exact_span_scores():
    context, names, database, species, _ = sequence_fixture()
    result = validate_local_profiles(context, names, database, species)
    assert result["accepted"] and not result["orthology_inferred"]
    assert result["parameters"]["database_residues"] == 2200
    assert result["parameters"]["database_genes"] == 7
    assert result["scorer"].startswith("original_cpu")
    assert len(result["profiles"]) == 2 and len(result["hits"]) == 6
    assert {row["target_gene"] for row in result["hits"]} == set(names) - {"outside"}
    from orthohmm.search.matrices import get_ka_params
    lam, k = get_ka_params("BLOSUM62")
    for row in result["hits"]:
        profile = result["profiles"]["ab".index(row["profile_side"])]
        assert row["target_gene"] not in profile["members"]
        assert profile["sequence_weighting"] == "none"
        assert row["supported"] and row["span"]["score"] == row["raw_score"]
        assert row["evalue"] == pytest.approx(k * profile["length"] * 2200 * np.exp(-lam * row["raw_score"]))
        assert len(profile["match_emissions"]) == profile["length"]
    assert result["supported_species_a"] == result["supported_species_b"] == ["sp1", "sp2", "sp3"]
    json.dumps(result, allow_nan=False)


def test_source_order_and_global_table_order_preserve_semantic_evidence():
    context, names, database, species, sequences = sequence_fixture()
    first = validate_local_profiles(context, names, database, species)
    other = deepcopy(context)
    other["genes_a"].reverse()
    other["genes_b"].reverse()
    names.reverse()
    reversed_db = species_sequences(names, [sequences[gene] for gene in names], "complete")
    assert validate_local_profiles(other, names, reversed_db, species) == first


def test_unrelated_sequences_do_not_nominate_a_joint_tree():
    context, names, database, species, _ = sequence_fixture(unrelated=True)
    result = validate_local_profiles(context, names, database, species)
    assert not result["accepted"] and result["reasons"]
    assert not any(row["supported"] for row in result["hits"])


def test_legitimate_asymmetric_fragments_are_not_automatically_rejected():
    context, names, _, species, _ = sequence_fixture()
    values = ["A" * (100 if gene.startswith("a") else 20) for gene in names]
    database = species_sequences(names, values, "complete")
    result = validate_local_profiles(context, names, database, species)
    assert result["accepted"]
    assert all(row["span"]["category"] == "asymmetric_coverage" for row in result["hits"])


def test_significant_partial_shared_regions_do_not_pass_the_alignment_gate():
    context, names, _, species, _ = sequence_fixture()
    values = ["A" * 200 if gene.startswith("a") else "C" * 80 + "A" * 60 + "C" * 80 for gene in names]
    database = species_sequences(names, values, "complete")
    result = validate_local_profiles(context, names, database, species)
    assert not result["accepted"]
    assert all(row["significant"] and not row["supported"] for row in result["hits"])
    assert all(row["reason"] == "aligned_coverage" for row in result["hits"])
    assert all(row["span"]["category"] == "partial_shared_region" for row in result["hits"])


def test_two_supported_genes_from_one_species_do_not_satisfy_profile_gate():
    context, names, _, species, sequences = sequence_fixture()
    species["b2"] = "sp1"
    sequences["b3"] = "W" * 200
    database = species_sequences(names, [sequences[gene] for gene in names], "complete")
    result = validate_local_profiles(context, names, database, species)
    assert result["supported_genes_b"] == ["b1", "b2"]
    assert result["supported_species_b"] == ["sp1"]
    assert not result["accepted"] and "insufficient_profile_species_b" in result["reasons"]


def test_uninformative_source_profile_is_an_explicit_abstention():
    context, names, _, species, _ = sequence_fixture()
    database = species_sequences(names, ["X" * 200 for gene in names], "complete")
    result = validate_local_profiles(context, names, database, species)
    assert result["reasons"] == ["missing_informative_profile"]
    assert not result["accepted"] and not result["hits"]


def test_original_scalar_fallback_is_explicit_and_score_equivalent(monkeypatch):
    from orthohmm.search import viterbi
    args = sequence_fixture()[:4]
    expected = validate_local_profiles(*args)
    def unavailable(*args, **kwargs):
        raise RuntimeError("multipair backend unavailable")
    monkeypatch.setattr(viterbi, "batch_viterbi_multipair_c", unavailable)
    actual = validate_local_profiles(*args)
    assert actual.pop("scorer") == "original_cpu_scalar_fallback"
    expected.pop("scorer")
    assert actual == expected


def test_no_silent_backend_replacement_when_original_cpu_library_missing(monkeypatch):
    from orthohmm.search import viterbi
    monkeypatch.setattr(viterbi, "is_c_available", lambda: False)
    with pytest.raises(RuntimeError, match="original CPU"):
        validate_local_profiles(*sequence_fixture()[:4])


@pytest.mark.parametrize("kind", ["duplicate", "missing", "overlap", "too_small", "one_species", "table_order", "offsets"])
def test_invalid_context_or_database_rejected_before_scoring(kind):
    context, names, database, species, _ = sequence_fixture()
    if kind == "duplicate":
        context["genes_a"][1] = "a1"
    elif kind == "missing":
        context["genes_a"][1] = "missing"
        species["missing"] = "sp2"
    elif kind == "overlap":
        context["genes_b"][1] = "a1"
    elif kind == "too_small":
        context["genes_a"].pop()
    elif kind == "one_species":
        for gene in context["genes_a"]:
            species[gene] = "sp1"
    elif kind == "table_order":
        names = list(reversed(names))
    else:
        database.offsets[1] += 1
    with pytest.raises(ValueError):
        validate_local_profiles(context, names, database, species)


def workflow_arguments(tmp_path, unrelated=False):
    context, names, database, species, sequences = sequence_fixture(unrelated=unrelated)
    path = tmp_path / "baseline_species.nwk"
    path.write_text("[&R] ((sp1,sp2),sp3);\n")
    config = PhylogenyConfig(mode="reconcile", species_tree_mode="infer", species_tree_rooting="min_variance",
        root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy",
        aligner=shutil.which("mafft") or "mafft", tree_builder=shutil.which("FastTree") or "FastTree")
    original = [(f"{side}{i}", f"{side}{j}") for side in "ab" for i in range(1, 4) for j in range(i + 1, 4)]
    return [context, names, database, species, sequences, path, ["sp1.faa", "sp2.faa", "sp3.faa"],
            config, tmp_path / "joint", [], original]


def test_profile_rejection_skips_tree_and_preserves_receipt(tmp_path, monkeypatch):
    import orthohmm.phylogeny_pipeline as pipeline
    def forbidden(*args, **kwargs):
        raise AssertionError("Rejected profile must not infer a tree")
    monkeypatch.setattr(pipeline, "_reconcile_family", forbidden)
    args = workflow_arguments(tmp_path, unrelated=True)
    result = evaluate_local_context(*args)
    assert result["status"] == "rejected_by_profiles" and not result["tree_inferred"]
    assert json.loads((args[8] / "proposal.json").read_text()) == result
    with pytest.raises(FileExistsError):
        evaluate_local_context(*args)


@pytest.mark.parametrize("kind", ["sequence", "config", "constraint", "baseline_pair", "context_id"])
def test_invalid_joint_context_rejected_before_output_creation(tmp_path, kind):
    args = workflow_arguments(tmp_path)
    if kind == "sequence":
        args[4]["a1"] = "A" * 200
    elif kind == "config":
        args[7] = replace(args[7], gene_tree_support_mode="internal_splits")
    elif kind == "constraint":
        args[9] = [{"source_genes": ["a1"], "target_genes": ["b1"]}]
    elif kind == "baseline_pair":
        args[10].append(("a1", "outside"))
    else:
        args[0]["context_id"] = "../escape"
    with pytest.raises(ValueError):
        evaluate_local_context(*args)
    assert not args[8].exists()


@pytest.mark.parametrize("stage", ["profile", "tree"])
def test_failed_context_keeps_a_terminal_failure_receipt(tmp_path, monkeypatch, stage):
    import orthohmm.local_rescue as local
    import orthohmm.phylogeny_pipeline as pipeline
    def failed(*args, **kwargs):
        raise RuntimeError("controlled tool failure")
    if stage == "profile":
        monkeypatch.setattr(local, "validate_local_profiles", failed)
    else:
        monkeypatch.setattr(pipeline, "_reconcile_family", failed)
    args = workflow_arguments(tmp_path)
    with pytest.raises(RuntimeError, match="controlled tool failure"):
        evaluate_local_context(*args)
    receipt = json.loads((args[8] / "proposal.json").read_text())
    assert receipt["status"] == ("profile_validation_failed" if stage == "profile" else "joint_tree_failed")
    assert not receipt["accepted"] and not receipt["tree_inferred"]
    assert receipt["error"] == {"type": "RuntimeError", "message": "controlled tool failure"}
    assert not receipt["families_merged"]


@pytest.mark.parametrize("value,label", [(float("inf"), "positive_infinity"),
    (float("-inf"), "negative_infinity"), (float("nan"), "nan")])
def test_legacy_constraint_diagnostics_are_json_safe_without_changing_runtime_constraints(tmp_path, monkeypatch, value, label):
    import orthohmm.local_rescue as local
    args = workflow_arguments(tmp_path, unrelated=True)
    constraint = {"source_genes": ["a1"], "target_genes": ["a2"], "margin": value,
                  "nested_diagnostic": [value, {"value": value}], "finite_diagnostic": 1.5}
    args[9] = [constraint]
    # A profile rejection still has to serialize its original constraint audit.
    result = evaluate_local_context(*args)
    encoded = {"nonfinite_constraint_diagnostic": label}
    assert result["constraints"][0]["margin"] == encoded
    assert result["constraints"][0]["nested_diagnostic"] == [encoded, {"value": encoded}]
    assert result["constraints"][0]["finite_diagnostic"] == 1.5
    assert result["constraints"][0]["source_genes"] == ["a1"]
    assert math.isnan(constraint["margin"]) if label == "nan" else constraint["margin"] == value
    assert json.loads((args[8] / "proposal.json").read_text()) == result
    assert not result["tree_inferred"] and not result["accepted"]
    assert local._context_constraints(args[9], [args[0]["genes_a"], args[0]["genes_b"]])[0] is constraint


def test_scoring_nonfinite_is_not_hidden_by_constraint_diagnostic_encoding(tmp_path, monkeypatch):
    import orthohmm.local_rescue as local
    args = workflow_arguments(tmp_path, unrelated=True)
    monkeypatch.setattr(local, "validate_local_profiles", lambda *a, **kw: {"accepted": False, "score": float("nan")})
    with pytest.raises(ValueError, match="JSON compliant"):
        evaluate_local_context(*args)


@pytest.mark.skipif(not shutil.which("mafft") or not shutil.which("FastTree"), reason="Requires real tree tools")
def test_original_infinite_diagnostic_reaches_actual_reconciliation_unchanged(tmp_path, monkeypatch):
    import orthohmm.phylogeny_pipeline as pipeline
    args = workflow_arguments(tmp_path)
    original = pipeline.apply_satellite_membership_constraints
    constraint = {"source_genes": ["a1"], "target_genes": ["a2"], "margin": float("inf")}
    args[9] = [constraint]
    calls = []
    def observed(outcomes, constraints, *a, **kw):
        calls.append(constraints)
        return original(outcomes, constraints, *a, **kw)
    monkeypatch.setattr(pipeline, "apply_satellite_membership_constraints", observed)
    result = evaluate_local_context(*args)
    assert result["tree_inferred"] and result["status"] == "joint_tree_evaluated"
    assert calls == [[constraint]] and calls[0][0] is constraint
    assert calls[0][0]["margin"] == float("inf")
    assert result["constraints"][0]["margin"] == {"nonfinite_constraint_diagnostic": "positive_infinity"}
    json.dumps(result, allow_nan=False)


def proposal_fixture(context, accepted=True):
    return {"context": context, "accepted": accepted, "status": "joint_tree_evaluated",
            "tree_inferred": True, "profile_validation": {"accepted": accepted},
            "joint_admission": {"accepted": accepted}}


def candidate_fixture(tmp_path):
    context, names, _, _, _ = sequence_fixture()
    original = tmp_path / "original.txt"
    raw = b"b3 b2 b1\r\n\r\n  a1\ta2 a3  \r\noutside\r\n"
    original.write_bytes(raw)
    constraints = [{"source_genes": ["a1"], "target_genes": ["a2", "a3"]}]
    return context, names, original, tmp_path / "effective.txt", constraints, raw


def test_rejected_or_no_proposals_preserve_original_bytes_exactly(tmp_path):
    context, names, original, output, constraints, raw = candidate_fixture(tmp_path)
    audit = write_effective_candidates(original, names, [proposal_fixture(context, False)], output, constraints)
    assert original.read_bytes() == output.read_bytes() == raw
    assert audit["original_bytes_preserved"] and audit["effective_families"] == 3
    assert not audit["accepted_contexts"] and not audit["final_native_inference"]
    with pytest.raises(FileExistsError):
        write_effective_candidates(original, names, [], output, constraints)


def test_only_accepted_whole_contexts_change_candidates_not_constraints(tmp_path):
    context, names, original, output, constraints, raw = candidate_fixture(tmp_path)
    before = deepcopy(constraints)
    audit = write_effective_candidates(original, names, [proposal_fixture(context)], output, constraints)
    assert audit["original_families"] == 3 and audit["effective_families"] == 2
    assert audit["genes"] == 7 and audit["constraints_preserved"] == 1
    assert output.read_text() == "a1 a2 a3 b1 b2 b3\noutside\n"
    assert original.read_bytes() == raw and constraints == before
    assert not audit["final_native_inference"]


@pytest.mark.parametrize("kind", ["partial_family", "overlap", "unfinished", "constraint", "duplicate_ownership", "missing_gene"])
def test_invalid_candidate_rewrite_leaves_output_absent(tmp_path, kind):
    context, names, original, output, constraints, _ = candidate_fixture(tmp_path)
    proposals = [proposal_fixture(context)]
    if kind == "partial_family":
        context["genes_a"].pop()
    elif kind == "overlap":
        proposals.append(deepcopy(proposals[0]))
    elif kind == "unfinished":
        proposals[0]["tree_inferred"] = False
    elif kind == "constraint":
        constraints[0]["target_genes"] = ["b1"]
    elif kind == "duplicate_ownership":
        original.write_text(original.read_text() + "a1\n")
    else:
        names.pop()
    with pytest.raises(ValueError):
        write_effective_candidates(original, names, proposals, output, constraints)
    assert not output.exists()


@pytest.mark.skipif(not shutil.which("mafft") or not shutil.which("FastTree"), reason="Requires real tree tools")
def test_fresh_real_mafft_fasttree_joint_reconciliation_and_exact_cache_reuse(tmp_path, monkeypatch):
    from orthohmm.phylogeny import parse_species_tree
    from orthohmm.phylogeny_pipeline import _load_reusable_raw_trees, _reconcile_family
    args = workflow_arguments(tmp_path)
    result = evaluate_local_context(*args)
    assert result["tree_inferred"] and result["status"] == "joint_tree_evaluated"
    assert result["profile_validation"]["accepted"]
    assert not result["families_merged"] and not result["final_native_inference"]
    directory = args[8]
    for suffix in ("raw", "rooted", "reconciled"):
        assert (directory / "gene_trees" / f"LocalContext000001.{suffix}.nwk").stat().st_size
    checkpoint = json.loads((directory / "checkpoints/LocalContext000001.json").read_text())
    assert checkpoint["genes"] == sorted(args[0]["genes_a"] + args[0]["genes_b"])
    assert len(result["token_to_gene"]) == 6
    cache = _load_reusable_raw_trees(directory)
    validated = parse_species_tree(args[5], args[6])
    genes = tuple(checkpoint["genes"])
    reused = _reconcile_family("Family0000001", genes, args[4], args[3], validated.tree,
        result["species_tree"]["sha256"], args[7], tmp_path / "final", cache)
    assert reused.checkpoint_hit and reused.remapped_checkpoint_hit
    assert reused.reconciliation is not None
    assert set(reused.genes) == set(genes)
    assert json.loads((directory / "proposal.json").read_text()) == result
    import orthohmm.phylogeny_pipeline as pipeline
    commands = []
    original = pipeline._run_to_file
    def observed(command, path, role):
        commands.append(command)
        return original(command, path, role)
    monkeypatch.setattr(pipeline, "_run_to_file", observed)
    changed = dict(args[4])
    changed["a1"] = "M" + changed["a1"][1:]
    if changed["a1"] == args[4]["a1"]:
        changed["a1"] = "A" + changed["a1"][1:]
    fresh = _reconcile_family("Family0000002", genes, changed, args[3], validated.tree,
        result["species_tree"]["sha256"], args[7], tmp_path / "changed", cache)
    assert not fresh.checkpoint_hit and not fresh.remapped_checkpoint_hit
    assert len(commands) == 2
    raw = directory / "gene_trees/LocalContext000001.raw.nwk"
    raw.write_text(raw.read_text() + "\n")
    assert _load_reusable_raw_trees(directory) == {}


@pytest.mark.skipif(not shutil.which("mafft") or not shutil.which("FastTree"), reason="Requires real tree tools")
def test_admitted_context_runs_through_ordinary_final_inferred_species_pipeline(tmp_path):
    from benchmark_tools.verify_frozen_phylogeny_install import fixture
    from orthohmm.phylogeny_pipeline import _load_sequence_data, run_phylogeny_stage
    from orthohmm.search.profile_expansion import load_global_sequence_database
    fasta = tmp_path / "input"
    fasta.mkdir()
    for taxon, genes in fixture().items():
        (fasta / f"{taxon}.faa").write_text("".join(f">{gene}\n{seq}\n" for gene, seq in genes.items()))
    files = sorted(path.name for path in fasta.iterdir())
    sequences, species = _load_sequence_data(str(fasta), files)
    names = sorted(sequences)
    groups = [[gene for gene in names if gene.endswith(f"family{index}")] for index in range(3)]
    groups.append([gene for gene in names if gene.endswith("duplicate")])
    original = tmp_path / "original_candidates.txt"
    original.write_text("".join(" ".join(group) + "\n" for group in groups))
    baseline = tmp_path / "baseline"
    cluster_path = baseline / "orthohmm_working_res/orthohmm_edges_clustered.txt"
    cluster_path.parent.mkdir(parents=True)
    shutil.copyfile(original, cluster_path)
    config = PhylogenyConfig(mode="reconcile", species_tree_mode="infer", species_tree_rooting="min_variance",
        root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy",
        aligner=shutil.which("mafft"), tree_builder=shutil.which("FastTree"))
    control = run_phylogeny_stage(str(fasta), str(baseline), files, config, 1)
    assert control.species_tree_families > 0 and not control.species_tree_checkpoint_hit
    phylo = baseline / "orthohmm_phylogeny"
    with (phylo / "orthohmm_pairwise_orthologs_confidence.tsv").open() as stream:
        pairs = [(row["gene_a"], row["gene_b"]) for row in csv.DictReader(stream, delimiter="\t")
                 if row["confidence"] == "high"]
    context = {"context_id": "LocalContext000001", "family_a": 2, "family_b": 3,
               "genes_a": groups[2], "genes_b": groups[3]}
    joint = tmp_path / "joint"
    proposal = evaluate_local_context(context, names, load_global_sequence_database(str(fasta), files, names),
        species, sequences, phylo / "species_tree.rooted.nwk", files, config, joint, [], pairs)
    assert proposal["accepted"] and proposal["joint_admission"]["novel_cross_high_pairs"]
    effective = tmp_path / "effective_candidates.txt"
    audit = write_effective_candidates(original, names, [proposal], effective, [])
    assert audit["effective_families"] == 3 and len(audit["accepted_contexts"]) == 1
    final = tmp_path / "final"
    final_clusters = final / "orthohmm_working_res/orthohmm_edges_clustered.txt"
    final_clusters.parent.mkdir(parents=True)
    shutil.copyfile(effective, final_clusters)
    for subdir in ("gene_trees", "checkpoints"):
        shutil.copytree(joint / subdir, final / "orthohmm_phylogeny" / subdir)
    actual = run_phylogeny_stage(str(fasta), str(final), files, config, 1)
    assert actual.candidate_families == 3 and actual.reconciled_families == 1
    assert actual.checkpoint_hits == actual.remapped_checkpoint_hits == 1
    assert actual.species_tree_families > 0 and not actual.species_tree_checkpoint_hit
    root_groups = [line.split() for line in final_clusters.read_text().splitlines()]
    assigned = [gene for group in root_groups for gene in group]
    assert len(assigned) == len(set(assigned)) == len(names) and set(assigned) == set(names)
    root_owner = {gene: index for index, group in enumerate(root_groups) for gene in group}
    effective_owner = {gene: index for index, group in enumerate(effective.read_text().splitlines()) for gene in group.split()}
    with (final / "orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv").open() as stream:
        final_pairs = list(csv.DictReader(stream, delimiter="\t"))
    assert final_pairs and len(final_pairs) == actual.ortholog_pairs
    for row in final_pairs:
        a, b = row["gene_a"], row["gene_b"]
        assert species[a] != species[b] and root_owner[a] == root_owner[b]
        assert effective_owner[a] == effective_owner[b]

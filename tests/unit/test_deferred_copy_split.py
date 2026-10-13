from argparse import Namespace
import json
from types import SimpleNamespace

import numpy as np
import pytest

from orthohmm import args_processing, refinement
from orthohmm.helpers import StartStep, SubstitutionMatrix
from orthohmm.orthohmm import execute
from orthohmm.parser import create_parser


def refine(size, maximum=0, dataset_species=78, scope="both"):
    species = np.array([0] * 14 + [1 + i % (dataset_species - 1) for i in range(size - 14)]
                       + list(range(dataset_species)), dtype=np.int32)
    clusters = [list(range(size))] + [[i] for i in range(size, len(species))]
    empty = np.array([], dtype=np.int32)
    return refinement.refine_cluster_indices(clusters, empty, empty, np.array([]),
        np.array([0], dtype=np.int32), np.array([14], dtype=np.int32), species,
        rbnh_scores=np.array([2.]), defer_copy_split_max_size=maximum, defer_copy_split_scope=scope)


@pytest.mark.parametrize("size", [70, 149, 150, 500])
def test_bounded_copy_rich_families_are_preserved_for_phylogeny(size):
    baseline = refine(size)
    deferred = refine(size, 500)
    assert list(range(size)) not in baseline
    assert list(range(size)) in deferred
    assert {g for group in baseline for g in group} == {g for group in deferred for g in group}
    assert len([g for group in deferred for g in group]) == size + 78
    assert deferred == refine(size, np.int32(500))


def test_larger_families_keep_baseline_splitting_and_explicit_zero_is_default():
    assert refine(501, 500) == refine(501)
    assert refine(150, 0) == refine(150)


def test_narrow_panel_component_splitting_is_also_deferred_within_bound():
    baseline = refine(70, dataset_species=12)
    deferred = refine(70, 500, dataset_species=12)
    assert list(range(70)) not in baseline
    assert list(range(70)) in deferred


def test_deferral_does_not_disable_noncopy_degree_splitting():
    cluster = list(range(40))
    species = np.array([0] * 20 + [1] * 20, dtype=np.int32)
    empty = np.array([], dtype=np.int32)
    q, t = np.zeros(5, dtype=np.int32), np.arange(1, 6, dtype=np.int32)
    args = ([cluster], empty, empty, np.array([]), q, t, species)
    baseline = refinement.refine_cluster_indices(*args)
    assert baseline != [cluster]
    assert refinement.refine_cluster_indices(*args, defer_copy_split_max_size=500) == baseline


@pytest.mark.parametrize("value", [-1, True, np.bool_(False), 1.5, "500", None])
def test_invalid_deferral_sizes_are_not_coerced(value):
    empty = np.array([], dtype=np.int32)
    with pytest.raises(ValueError, match="nonnegative integer"):
        refinement.refine_cluster_indices([], empty, empty, np.array([]), empty, empty, empty,
                                          defer_copy_split_max_size=value)


def test_named_preset_does_not_mutate_defaults_or_existing_qfo_preset():
    assert refinement.resolve_refinement_profile("default") == {}
    qfo = refinement.resolve_refinement_profile("qfo")
    preset = refinement.resolve_refinement_profile("phylogeny_deferred")
    assert preset == {"defer_copy_split_max_size": 500}
    preset["defer_copy_split_max_size"] = 999
    assert refinement.resolve_refinement_profile("phylogeny_deferred") == {"defer_copy_split_max_size": 500}
    assert refinement.resolve_refinement_profile("qfo") == qfo


def cli_args(tmp_path):
    (tmp_path / "A.faa").write_text(">a\nACDE\n")
    (tmp_path / "B.faa").write_text(">b\nACDE\n")
    return Namespace(fasta_directory=str(tmp_path), output_directory=str(tmp_path), phmmer=None,
        cpu=1, single_copy_threshold=.5, mcl="mcl", inflation_value=1.5, start=None, stop=None,
        substitution_matrix="BLOSUM62", evalue=1e-4, search_mode="builtin", clustering="leiden",
        cpm_resolution=.1, refinement_profile="phylogeny_deferred", accuracy_profile="high_sensitivity",
        phylogeny="reconcile", phylogeny_candidates="satellite_v2", phylogeny_root_rule="species_overlap",
        phylogeny_pair_rule="positive_paralogy", species_tree_mode="infer", species_tree_rooting="min_variance",
        aligner="/bin/true", tree_builder="/bin/true")


def test_native_cli_accepts_explicit_deferred_profile_in_required_context(tmp_path):
    args = cli_args(tmp_path)
    parsed = create_parser().parse_args([str(tmp_path), "--refinement_profile", "phylogeny_deferred"])
    assert parsed.refinement_profile == "phylogeny_deferred"
    result = args_processing.process_args(args)
    assert result["refinement_profile"] == "phylogeny_deferred"
    assert result["phylogeny_pair_rule"] == "positive_paralogy"


@pytest.mark.parametrize("field,value", [
    ("accuracy_profile", "standard"), ("phylogeny", "off"), ("phylogeny_candidates", "seed"),
    ("phylogeny_root_rule", "supported_children"), ("phylogeny_pair_rule", "lca")])
def test_cli_rejects_deferred_profile_without_required_phylogenetic_context(tmp_path, field, value):
    args = cli_args(tmp_path)
    setattr(args, field, value)
    with pytest.raises(SystemExit):
        args_processing.process_args(args)


def api_kwargs():
    return dict(accuracy_profile="high_sensitivity", refinement_profile="phylogeny_deferred",
        phylogeny="reconcile", phylogeny_candidates="satellite_v2", phylogeny_root_rule="species_overlap",
        phylogeny_pair_rule="positive_paralogy")


@pytest.mark.parametrize("field,value", [
    ("accuracy_profile", "standard"), ("phylogeny", "off"), ("phylogeny_candidates", "seed"),
    ("phylogeny_root_rule", "supported_children"), ("phylogeny_pair_rule", "lca"),
    ("start", StartStep.search_res)])
def test_programmatic_pipeline_enforces_context_before_search(tmp_path, field, value):
    kwargs = api_kwargs()
    start = value if field == "start" else None
    if field != "start":
        kwargs[field] = value
    with pytest.raises(ValueError, match="phylogeny_deferred refinement requires"):
        execute(str(tmp_path), str(tmp_path / "out"), None, 1, .5, "mcl", 1.5,
                start, None, SubstitutionMatrix.blosum62, 1e-4, **kwargs)


def test_pipeline_records_experimental_policy_and_accepts_case_normalization(tmp_path, monkeypatch):
    from orthohmm.search import engine
    class SearchReached(Exception):
        pass
    def reached(*args, **kwargs):
        raise SearchReached()
    monkeypatch.setattr(engine, "execute_builtin_search", reached)
    (tmp_path / "species.faa").write_text(">gene\nACDE\n")
    metrics = tmp_path / "metrics.json"
    kwargs = api_kwargs()
    kwargs["refinement_profile"] = "PHYLOGENY_DEFERRED"
    with pytest.raises(SearchReached):
        execute(str(tmp_path), str(tmp_path / "out"), None, 1, .5, "mcl", 1.5,
                None, None, SubstitutionMatrix.blosum62, 1e-4, metrics_json=str(metrics), **kwargs)
    metadata = json.loads(metrics.read_text())["metadata"]
    assert metadata["refinement_profile"] == "phylogeny_deferred"
    assert metadata["defer_copy_split_max_size"] == 500


def test_benchmark_wrapper_forwards_deferred_profile(tmp_path, monkeypatch):
    from benchmark_tools import benchmark_production
    fasta = tmp_path / "fasta"
    fasta.mkdir()
    (fasta / "A.faa").write_text(">a\nACDE\n")
    captured = {}
    def fake_run(command, **kwargs):
        captured["command"] = command
        with open(command[command.index("--metrics_json") + 1], "w") as stream:
            json.dump({"schema_version": 1, "status": "complete"}, stream)
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(benchmark_production.subprocess, "run", fake_run)
    monkeypatch.setattr(benchmark_production, "git_value", lambda *args: "")
    assert benchmark_production.main([str(fasta), str(tmp_path / "out"), str(tmp_path / "result.json"),
        "--cpu", "1", "--refinement-profile", "phylogeny_deferred"]) == 0
    command = captured["command"]
    assert command[command.index("--refinement_profile") + 1] == "phylogeny_deferred"


@pytest.mark.parametrize("size,branch", [(70, "medium"), (149, "medium"), (150, "large"), (500, "large")])
def test_scoped_deferral_keeps_only_the_prespecified_branch(size, branch):
    other = "large" if branch == "medium" else "medium"
    assert refine(size, 500, scope=branch) == refine(size, 500)
    assert refine(size, 500, scope=other) == refine(size)
    assert refine(size, scope=branch) == refine(size)


@pytest.mark.parametrize("scope", ["both", "large", "medium"])
def test_all_scopes_preserve_genes_order_determinism_and_larger_family_safeguard(scope):
    for size in (70, 149, 150, 500, 501):
        groups = refine(size, 500, scope=scope)
        assert groups == refine(size, 500, scope=scope)
        genes = [gene for group in groups for gene in group]
        assert len(genes) == len(set(genes)) == size + 78
        assert set(genes) == set(range(size + 78))
    assert refine(501, 500, scope=scope) == refine(501)


@pytest.mark.parametrize("dataset_species", [12, 49, 50, 78])
def test_unchanged_panel_eligibility_not_dataset_name_controls_large_deferral(dataset_species):
    baseline = refine(150, dataset_species=dataset_species)
    large = refine(150, 500, dataset_species=dataset_species, scope="large")
    if dataset_species < 50:
        assert large == baseline
    else:
        assert list(range(150)) not in baseline and list(range(150)) in large


def test_narrow_panel_large_scope_is_inactive_and_medium_reproduces_both():
    assert refine(70, 500, dataset_species=12, scope="large") == refine(70, dataset_species=12)
    assert refine(70, 500, dataset_species=12, scope="medium") == refine(70, 500, dataset_species=12)


@pytest.mark.parametrize("scope", [None, True, 0, "Large", "all", "", [], np.array(["large"])])
def test_invalid_scope_is_rejected_even_for_empty_input(scope):
    empty = np.array([], dtype=np.int32)
    with pytest.raises(ValueError, match="scope"):
        refinement.refine_cluster_indices([], empty, empty, np.array([]), empty, empty, empty,
                                          defer_copy_split_scope=scope)


@pytest.mark.parametrize("branch", ["large", "medium"])
def test_scoped_preset_is_explicit_and_does_not_mutate_existing_profiles(branch):
    profile = "phylogeny_deferred_" + branch
    assert refinement.resolve_refinement_profile(profile) == {
        "defer_copy_split_max_size": 500, "defer_copy_split_scope": branch}
    assert refinement.resolve_refinement_profile("default") == {}
    assert refinement.resolve_refinement_profile("phylogeny_deferred") == {"defer_copy_split_max_size": 500}


@pytest.mark.parametrize("branch", ["large", "medium"])
def test_native_cli_accepts_scoped_profile_only_in_required_context(tmp_path, branch):
    profile = "phylogeny_deferred_" + branch
    args = cli_args(tmp_path)
    args.refinement_profile = profile
    parsed = create_parser().parse_args([str(tmp_path), "--refinement_profile", profile])
    assert parsed.refinement_profile == profile
    assert args_processing.process_args(args)["refinement_profile"] == profile
    args.phylogeny_pair_rule = "lca"
    with pytest.raises(SystemExit):
        args_processing.process_args(args)


@pytest.mark.parametrize("branch", ["large", "medium"])
@pytest.mark.parametrize("field,value", [
    ("accuracy_profile", "standard"), ("phylogeny", "off"), ("phylogeny_candidates", "seed"),
    ("phylogeny_root_rule", "supported_children"), ("phylogeny_pair_rule", "lca"),
    ("start", StartStep.search_res)])
def test_scoped_programmatic_profiles_cannot_bypass_phylogenetic_guards(tmp_path, branch, field, value):
    kwargs = api_kwargs()
    kwargs["refinement_profile"] = "phylogeny_deferred_" + branch
    start = value if field == "start" else None
    if field != "start":
        kwargs[field] = value
    with pytest.raises(ValueError, match="phylogeny_deferred refinement requires"):
        execute(str(tmp_path), str(tmp_path / "out"), None, 1, .5, "mcl", 1.5,
                start, None, SubstitutionMatrix.blosum62, 1e-4, **kwargs)


@pytest.mark.parametrize("branch", ["large", "medium"])
def test_scoped_pipeline_records_actual_experimental_policy_before_search(tmp_path, monkeypatch, branch):
    from orthohmm.search import engine
    class SearchReached(Exception):
        pass
    def reached(*args, **kwargs):
        raise SearchReached()
    monkeypatch.setattr(engine, "execute_builtin_search", reached)
    (tmp_path / "species.faa").write_text(">gene\nACDE\n")
    metrics = tmp_path / "metrics.json"
    kwargs = api_kwargs()
    kwargs["refinement_profile"] = ("phylogeny_deferred_" + branch).upper()
    with pytest.raises(SearchReached):
        execute(str(tmp_path), str(tmp_path / "out"), None, 1, .5, "mcl", 1.5,
                None, None, SubstitutionMatrix.blosum62, 1e-4, metrics_json=str(metrics), **kwargs)
    metadata = json.loads(metrics.read_text())["metadata"]
    assert metadata["refinement_profile"] == "phylogeny_deferred_" + branch
    assert metadata["defer_copy_split_scope"] == branch and metadata["defer_copy_split_max_size"] == 500


@pytest.mark.parametrize("dataset_species", [12, 78])
@pytest.mark.parametrize("scope", ["large", "medium"])
def test_scope_is_forwarded_after_reciprocal_merge_paths(monkeypatch, dataset_species, scope):
    species = np.array([0, 1, 0, 1, *range(dataset_species)], dtype=np.int32)
    clusters = [[0, 1], [2, 3], *[[i] for i in range(4, len(species))]]
    query, target, weight = np.array([0]), np.array([2]), np.array([2.])
    observed = []
    def split(groups, *args, **kwargs):
        observed.append(kwargs)
        return groups
    monkeypatch.setattr(refinement, "_split_high_duplication_index_clusters", split)
    refinement.refine_cluster_indices(clusters, query, target, weight, query, target, species,
        rbnh_scores=weight, defer_copy_split_max_size=500, defer_copy_split_scope=scope)
    assert len(observed) == 1 and observed[0]["defer_copy_split_scope"] == scope
    assert observed[0]["defer_copy_split_max_size"] == 500
    assert observed[0].get("broad_mode", False) is (dataset_species >= 50)


@pytest.mark.parametrize("branch", ["large", "medium"])
def test_benchmark_wrapper_forwards_scoped_profile(tmp_path, monkeypatch, branch):
    from benchmark_tools import benchmark_production
    fasta = tmp_path / "fasta"
    fasta.mkdir()
    (fasta / "A.faa").write_text(">a\nACDE\n")
    captured = {}
    def fake_run(command, **kwargs):
        captured["command"] = command
        with open(command[command.index("--metrics_json") + 1], "w") as stream:
            json.dump({"schema_version": 1, "status": "complete"}, stream)
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(benchmark_production.subprocess, "run", fake_run)
    monkeypatch.setattr(benchmark_production, "git_value", lambda *args: "")
    profile = "phylogeny_deferred_" + branch
    assert benchmark_production.main([str(fasta), str(tmp_path / "out"), str(tmp_path / "result.json"),
        "--cpu", "1", "--refinement-profile", profile]) == 0
    command = captured["command"]
    assert command[command.index("--refinement_profile") + 1] == profile

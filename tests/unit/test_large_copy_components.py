from argparse import Namespace
import json
from types import SimpleNamespace

import numpy as np
import pytest

from orthohmm import args_processing, refinement
from orthohmm.helpers import StartStep, SubstitutionMatrix
from orthohmm.orthohmm import execute
from orthohmm.parser import create_parser


def inputs(size=150, species_count=78):
    species = np.array([0] * 14 + [1 + i % (species_count - 1) for i in range(size - 14)]
                       + list(range(species_count)), dtype=np.int32)
    clusters = [list(range(size)), *[[i] for i in range(size, len(species))]]
    return clusters, species


def refine(size=150, edges=(), scores=None, maximum=500, species_count=78, **kwargs):
    clusters, species = inputs(size, species_count)
    query = np.array([a for a, b in edges], dtype=np.int32)
    target = np.array([b for a, b in edges], dtype=np.int32)
    empty = np.array([], dtype=np.int32)
    return refinement.refine_cluster_indices(clusters, empty, empty, [], query, target, species,
        rbnh_scores=np.full(len(edges), 2.) if scores is None else scores,
        retain_large_copy_components_max_size=maximum, **kwargs)


def test_strong_components_are_candidates_and_weak_bridges_do_not_join_them():
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]
    groups = refine(edges=edges, scores=[2., 1.5, 1.49, 2., 2.])
    assert groups[:2] == [[0, 1, 2], [3, 4, 5]]
    assert groups[2:150 - 4] == [[i] for i in range(6, 150)]
    genes = [g for c in groups for g in c]
    assert sorted(genes) == list(range(228)) and len(genes) == len(set(genes))
    assert all(len(c) == 1 for c in refine(edges=edges, maximum=0))


@pytest.mark.parametrize("size", [150, 500, 501, 600])
def test_oversized_families_can_retain_only_bounded_components(size):
    assert refine(size, [(0, 1), (1, 2)])[0] == [0, 1, 2]
    assert sorted(g for c in refine(size, [(0, 1), (1, 2)]) for g in c) == list(range(size + 78))


@pytest.mark.parametrize("component_size", [499, 500, 501])
def test_component_bound_is_inclusive_and_larger_components_singletonize(component_size):
    edges = [(i, i + 1) for i in range(component_size - 1)]
    groups = refine(600, edges)
    if component_size <= 500:
        assert groups[0] == list(range(component_size))
    else:
        assert all(len(c) == 1 for c in groups)
    assert sorted(g for c in groups for g in c) == list(range(678))


def test_empty_graph_has_no_unsupported_whole_family_fallback():
    assert refine(600) == refine(600, maximum=0)


def test_cross_family_edges_and_self_loops_do_not_leak_members():
    groups = refine(edges=[(0, 1), (1, 150), (2, 2)])
    assert groups[0] == [0, 1] and [150] in groups and [2] in groups
    assert not any(150 in c and 1 in c for c in groups)


def test_components_are_canonical_under_edge_orientation_order_and_cluster_order():
    edges = [(4, 3), (0, 2), (1, 0), (4, 5), (0, 1)]
    expected = refine(edges=edges)
    assert expected == refine(edges=list(reversed(edges)))
    clusters, species = inputs()
    clusters[0].reverse()
    empty = np.array([], dtype=np.int32)
    groups = refinement.refine_cluster_indices(clusters, empty, empty, [],
        np.array([a for a, b in edges]), np.array([b for a, b in edges]), species,
        rbnh_scores=np.full(len(edges), 2.), retain_large_copy_components_max_size=500)
    assert groups == expected


@pytest.mark.parametrize("size,species_count", [(70, 12), (149, 78), (150, 49)])
def test_medium_and_narrow_panel_behavior_is_unchanged(size, species_count):
    assert refine(size, [(0, 14)], species_count=species_count) == refine(
        size, [(0, 14)], maximum=0, species_count=species_count)


def test_trace_distinguishes_retained_oversized_and_isolated_components():
    trace = []
    edges = [(i, i + 1) for i in range(500)] + [(501, 502)]
    groups = refine(600, edges, copy_split_trace=trace)
    assert len(trace) == 1 and trace[0]["deferred"] is False
    assert trace[0]["component_sizes"] == [501, 2, *[1] * 97]
    assert trace[0]["retained_component_sizes"] == [2]
    assert trace[0]["oversized_component_sizes"] == [501]
    assert trace[0]["isolated_genes"] == 97
    assert trace[0]["component_min_edge_weight"] == 1.5
    assert [501, 502] in groups
    baseline = []
    refine(600, edges, maximum=0, copy_split_trace=baseline)
    assert set(baseline[0]) == {"members", "large_copy_eligible", "medium_copy_eligible", "deferred",
                               "size", "maximum_species_copies", "dataset_species_count"}


@pytest.mark.parametrize("value", [-1, 501, True, np.bool_(False), 1.5, "500", None])
def test_component_bound_is_validated_before_empty_return(value):
    with pytest.raises(ValueError, match="nonnegative integer at most 500"):
        refinement.refine_cluster_indices([], [], [], [], [], [], [],
            retain_large_copy_components_max_size=value)


def test_component_retention_cannot_be_combined_with_whole_family_deferral():
    with pytest.raises(ValueError, match="cannot be combined"):
        refine(defer_copy_split_max_size=500)


@pytest.mark.parametrize("queries,targets,scores,message", [
    ([0], [], [2.], "aligned"), ([0], [1], [], "aligned"),
    ([0], [1], None, "aligned"), ([[0]], [[1]], [2.], "one-dimensional"),
    ([0], [1], [[2.]], "one-dimensional"),
    ([0.5], [1], [2.], "integer"), ([True], [1], [2.], "integer"),
    ([-1], [1], [2.], "outside"), ([0], [228], [2.], "outside"),
    ([0], [1], [float("nan")], "finite positive"),
    ([0], [1], [float("inf")], "finite positive"),
    ([0], [1], [0.], "finite positive"), ([0], [1], [-2.], "finite positive"),
])
def test_component_scoring_evidence_is_not_silently_truncated_or_coerced(queries, targets, scores, message):
    clusters, species = inputs()
    with pytest.raises(ValueError, match=message):
        refinement.refine_cluster_indices(clusters, [], [], [], queries, targets, species,
            rbnh_scores=scores, retain_large_copy_components_max_size=500)


@pytest.mark.parametrize("threshold", [0., -1., float("nan"), float("inf")])
def test_invalid_component_threshold_is_rejected(threshold):
    with pytest.raises(ValueError, match="finite positive"):
        refine(panel_copy_component_min_edge_weight=threshold)


@pytest.mark.parametrize("species_count", [12, 78])
@pytest.mark.parametrize("between_clusters", [False, True])
def test_component_policy_is_forwarded_through_all_merge_paths(monkeypatch, species_count, between_clusters):
    clusters, species = inputs(species_count=species_count)
    q, t, w = np.array([0]), np.array([150 if between_clusters else 1]), np.array([2.])
    captured = []
    monkeypatch.setattr(refinement, "_split_high_duplication_index_clusters",
                        lambda groups, *args, **kwargs: captured.append(kwargs) or groups)
    refinement.refine_cluster_indices(clusters, q, t, w, q, t, species,
        rbnh_scores=w, retain_large_copy_components_max_size=np.int32(500))
    assert len(captured) == 1 and captured[0]["retain_large_copy_components_max_size"] == 500


def api_kwargs():
    return dict(accuracy_profile="high_sensitivity", refinement_profile="phylogeny_components",
        phylogeny="reconcile", phylogeny_candidates="satellite_v2", phylogeny_root_rule="species_overlap",
        phylogeny_pair_rule="positive_paralogy", species_tree_mode="infer", species_tree_rooting="min_variance")


def cli_args(tmp_path):
    for name in ("A", "B"):
        (tmp_path / (name + ".faa")).write_text(">" + name + "\nACDE\n")
    return Namespace(fasta_directory=str(tmp_path), output_directory=str(tmp_path), phmmer=None,
        cpu=1, single_copy_threshold=.5, mcl="mcl", inflation_value=1.5, start=None, stop=None,
        substitution_matrix="BLOSUM62", evalue=1e-4, search_mode="builtin", clustering="leiden",
        cpm_resolution=.1, aligner="/bin/true", tree_builder="/bin/true", **api_kwargs())


@pytest.mark.parametrize("field,value", [
    ("accuracy_profile", "standard"), ("phylogeny", "off"), ("phylogeny_candidates", "seed"),
    ("phylogeny_root_rule", "supported_children"), ("phylogeny_pair_rule", "lca"),
    ("start", StartStep.search_res), ("satellite_membership_rule", "path_species"),
    ("species_tree_mode", "provided"), ("species_tree_rooting", "midpoint"),
])
def test_api_rejects_component_profile_without_required_native_context(tmp_path, field, value):
    kwargs = api_kwargs()
    start = value if field == "start" else None
    if field != "start":
        kwargs[field] = value
    with pytest.raises(ValueError, match="phylogeny_components refinement requires"):
        execute(str(tmp_path), str(tmp_path / "out"), None, 1, .5, "mcl", 1.5,
                start, None, SubstitutionMatrix.blosum62, 1e-4, **kwargs)


@pytest.mark.parametrize("field,value", [
    ("accuracy_profile", "standard"), ("phylogeny", "off"), ("phylogeny_candidates", "seed"),
    ("phylogeny_root_rule", "supported_children"), ("phylogeny_pair_rule", "lca"),
    ("start", "search_res"), ("satellite_membership_rule", "path_species"),
    ("species_tree_mode", "provided"), ("species_tree_rooting", "midpoint"),
])
def test_cli_enforces_component_context(tmp_path, field, value):
    args = cli_args(tmp_path)
    setattr(args, field, value)
    with pytest.raises(SystemExit):
        args_processing.process_args(args)


def test_component_profile_is_opt_in_and_presets_are_not_mutated(tmp_path):
    expected = {"retain_large_copy_components_max_size": 500}
    preset = refinement.resolve_refinement_profile("phylogeny_components")
    assert preset == expected
    preset.clear()
    assert refinement.resolve_refinement_profile("PHYLOGENY_COMPONENTS") == expected
    assert refinement.resolve_refinement_profile("default") == {}
    assert refinement.resolve_refinement_profile("phylogeny_deferred") == {"defer_copy_split_max_size": 500}
    args = cli_args(tmp_path)
    assert args_processing.process_args(args)["refinement_profile"] == "phylogeny_components"
    parsed = create_parser().parse_args([str(tmp_path), "--refinement_profile", "phylogeny_components"])
    assert parsed.refinement_profile == "phylogeny_components"


def test_pipeline_records_actual_component_policy_before_search(tmp_path, monkeypatch):
    from orthohmm.search import engine
    class SearchReached(Exception):
        pass
    def reached(*args, **kwargs):
        raise SearchReached()
    monkeypatch.setattr(engine, "execute_builtin_search", reached)
    (tmp_path / "species.faa").write_text(">gene\nACDE\n")
    kwargs = api_kwargs()
    kwargs["refinement_profile"] = "PHYLOGENY_COMPONENTS"
    metrics = tmp_path / "metrics.json"
    with pytest.raises(SearchReached):
        execute(str(tmp_path), str(tmp_path / "out"), None, 1, .5, "mcl", 1.5,
                None, None, SubstitutionMatrix.blosum62, 1e-4, metrics_json=str(metrics), **kwargs)
    metadata = json.loads(metrics.read_text())["metadata"]
    assert metadata["refinement_profile"] == "phylogeny_components"
    assert metadata["retain_large_copy_components_max_size"] == 500
    assert metadata["large_copy_component_min_edge_weight"] == 1.5
    assert metadata["oversized_large_copy_component_policy"] == "singletons"
    assert "defer_copy_split_max_size" not in metadata


def test_benchmark_wrapper_forwards_component_profile(tmp_path, monkeypatch):
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
        "--cpu", "1", "--refinement-profile", "phylogeny_components"]) == 0
    command = captured["command"]
    assert command[command.index("--refinement_profile") + 1] == "phylogeny_components"


@pytest.mark.parametrize("numeric", [False, True])
def test_native_cluster_file_boundary_uses_actual_scored_edges_and_preserves_all_genes(tmp_path, numeric):
    from orthohmm.helpers import IndexedEdges
    from orthohmm.orthohmm import _refine_cluster_file
    clusters, species = inputs()
    names = ["gene" + str(i) for i in range(len(species))]
    lengths = [{"name": name, "spp": str(s)} for name, s in zip(names, species)]
    working = tmp_path / "orthohmm_working_res"
    working.mkdir()
    path = working / "orthohmm_edges_clustered.txt"
    path.write_text("".join(" ".join(names[i] for i in cluster) + "\n" for cluster in clusters))
    edges = (IndexedEdges(names, np.array([0, 1, 2]), np.array([1, 2, 3]), np.array([2., 1.5, 1.49]))
             if numeric else {frozenset((names[a], names[b])): w
                              for a, b, w in ((0, 1, 2.), (1, 2, 1.5), (2, 3, 1.49))})
    _refine_cluster_file(str(tmp_path), lengths, {}, edges, 1e-4, refinement_profile="phylogeny_components")
    groups = [line.split() for line in path.read_text().splitlines()]
    assert groups[0] == ["gene0", "gene1", "gene2"] and ["gene3"] in groups
    assert sorted(g for c in groups for g in c) == sorted(names)
    assert len([g for c in groups for g in c]) == len(names)


@pytest.mark.parametrize("field,value", [("search_mode", "phmmer"), ("clustering", "mcl")])
def test_component_native_pipeline_cannot_bypass_hmm_and_leiden_requirements(tmp_path, field, value):
    kwargs = api_kwargs()
    kwargs[field] = value
    with pytest.raises(ValueError, match="high_sensitivity accuracy requires"):
        execute(str(tmp_path), str(tmp_path / "out"), None, 1, .5, "mcl", 1.5,
                None, None, SubstitutionMatrix.blosum62, 1e-4, **kwargs)

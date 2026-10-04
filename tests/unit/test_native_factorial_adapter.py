import ast
from dataclasses import asdict
import itertools
from pathlib import Path

import pytest

from benchmark_tools import native_factorial_adapter as adapter
from orthohmm import orthohmm as pipeline


CELLS = tuple("p%d_c%d_r%d" % v for v in itertools.product(range(2), repeat=3))


@pytest.mark.parametrize("cell", CELLS)
def test_exact_ast_scope_for_every_factor_combination(cell):
    source = Path(pipeline.__file__).read_bytes()
    original = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == "_execute")
    observed = adapter.execution_tree(source, cell).body[0]
    changed = cell[4] == "1" and cell[7] == "0"
    if changed:
        gate = next(n for n in ast.walk(original) if isinstance(n, ast.If) and ast.dump(n.test) == ast.dump(adapter.GATE))
        gate.test = gate.test.values[1]
    assert ast.dump(observed) == ast.dump(original)


@pytest.mark.parametrize("cell", CELLS)
def test_profile_off_keeps_sensitive_search_graph_seed_and_production_globals(cell):
    original_execute = pipeline._execute
    original_resolver = pipeline.resolve_accuracy_profile
    old = asdict(original_resolver("high_sensitivity"))
    adapted, report = adapter.entrypoint(pipeline, cell)
    active = adapted.__globals__["resolve_accuracy_profile"]("high_sensitivity")
    assert asdict(active) == dict(old, profile_expansion=cell[1] == "1")
    assert (active.kmer_k, active.max_candidates_per_query, active.multipass_graph, active.leiden_seed) == (4, 100, True, 4)
    assert pipeline._execute is original_execute and pipeline.resolve_accuracy_profile is original_resolver
    assert asdict(original_resolver("high_sensitivity")) == old
    assert report["candidate_gate_decoupled"] == (cell[4] == "1" and cell[7] == "0")
    with pytest.raises(ValueError, match="substitute standard"):
        adapted.__globals__["resolve_accuracy_profile"]("standard")


@pytest.mark.parametrize("cell", [None, "p2_c1_r0", "p0c1r0", "p1_c0_r0\n", 123, "p1_c0_r0_extra"])
def test_invalid_factor_identity_is_rejected(cell):
    with pytest.raises(ValueError):
        adapter.factors(cell)


def test_frozen_source_change_refused_before_compilation():
    source = Path(pipeline.__file__).read_bytes()
    with pytest.raises(ValueError, match="exact frozen"):
        adapter.execution_tree(source + b"\n", "p0_c1_r0")


def test_factor_controls_are_independent_and_cannot_mutate_each_other():
    off, _ = adapter.entrypoint(pipeline, "p0_c0_r0")
    on, _ = adapter.entrypoint(pipeline, "p1_c1_r1")
    assert off.__globals__ is not on.__globals__
    assert off.__globals__["resolve_accuracy_profile"]("high_sensitivity").profile_expansion is False
    assert on.__globals__["resolve_accuracy_profile"]("high_sensitivity").profile_expansion is True


@pytest.mark.parametrize("change", [
    {"accuracy_profile": "standard"}, {"cpm_resolution": .2}, {"start": "search_res"},
    {"phylogeny": "reconcile"}, {"phylogeny_candidates": "seed"}, {"species_tree_mode": "supplied"},
    {"evalue_threshold": .001}, {"refinement_profile": "qfo"},
])
def test_changed_factorial_recipe_rejected_before_pipeline(change):
    entry, _ = adapter.entrypoint(pipeline, "p0_c1_r0")
    kwargs = dict(accuracy_profile="high_sensitivity", search_mode="builtin", clustering="leiden",
        cpm_resolution=.1, refinement_profile="default", evalue_threshold=.0001, start=None,
        phylogeny="off", phylogeny_candidates="satellite_v2", species_tree_mode="infer",
        species_tree_rooting="min_variance", phylogeny_root_rule="species_overlap",
        phylogeny_pair_rule="positive_paralogy")
    kwargs.update(change)
    with pytest.raises(ValueError, match="explicit frozen"):
        entry.__globals__["_execute"](**kwargs)

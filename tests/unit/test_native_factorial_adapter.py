import ast
import csv
from contextlib import contextmanager
from dataclasses import asdict
import itertools
import hashlib
import json
from pathlib import Path

import numpy as np
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


@pytest.mark.parametrize("cell", CELLS)
def test_candidate_gate_with_actual_nonzero_merges_independent_of_reconciliation(tmp_path, cell):
    working = tmp_path / "orthohmm_working_res"
    working.mkdir()
    names = ["g%d" % i for i in range(19)]
    clusters = [[i] for i in range(9)] + [list(range(9, 19))]
    seed = "".join(" ".join(names[i] for i in group) + "\n" for group in clusters)
    path = working / "orthohmm_edges_clustered.txt"
    path.write_text(seed)
    queries, targets = [], []
    for satellite in range(9):
        for anchor in range(9, 19):
            queries.extend((satellite, anchor))
            targets.extend((anchor, satellite))
    hits = (np.array(queries, dtype=np.int32), np.array(targets, dtype=np.int32), np.ones(180))
    class Metrics:
        def __init__(self):
            self.stages, self.metadata, self.counts = [], {}, {}
        @contextmanager
        def stage(self, label):
            self.stages.append(label)
            yield
        def add_metadata(self, **kwargs):
            self.metadata.update(kwargs)
        def add_counts(self, **kwargs):
            self.counts.update(kwargs)
    flags = adapter.factors(cell)
    tree = adapter.execution_tree(Path(pipeline.__file__).read_bytes(), cell)
    gate = next(n for n in ast.walk(tree) if isinstance(n, ast.If) and any(
        isinstance(child, ast.Call) and isinstance(child.func, ast.Name)
        and child.func.id == "_expand_phylogeny_candidates" for child in ast.walk(n)))
    metrics = Metrics()
    namespace = dict(pipeline.__dict__, metrics=metrics, output_directory=str(tmp_path),
        gene_names=names, gene_to_species=np.array([i % 3 for i in range(19)]), accuracy_hits=hits,
        phylogeny="reconcile" if flags["reconciliation"] else "off",
        phylogeny_candidates="satellite_v2" if flags["candidate_expansion"] else "seed")
    exec(compile(ast.fix_missing_locations(ast.Module(body=[gate], type_ignores=[])), "factorial-gate", "exec"), namespace)
    if flags["candidate_expansion"]:
        assert metrics.stages == ["phylogeny_candidates"]
        assert metrics.counts == {"phylogeny_seed_families": 10, "phylogeny_candidate_merges": 8}
        groups = [line.split() for line in path.read_text().splitlines()]
        assert sorted(map(len, groups)) == [1, 18]
        assert [group for group in groups if len(group) == 1] == [["g8"]]
        assert namespace["candidate_membership_constraints"] is not None
    else:
        assert metrics.stages == [] and metrics.counts == {} and path.read_text() == seed


def test_actual_eight_native_diagnostics_from_copied_raw_artifacts():
    root = Path(__file__).resolve().parents[2]
    copied = root / "benchmark_tools/results/native_factorial_adapter_diagnostic_20261004"
    report = json.loads((copied / "probe.json").read_text())
    assert report["status"] == "eight_native_factorial_diagnostics_complete"
    assert report["genes"] == 26 and report["hits"] == 116
    assert report["diagnostic_only"] is True and report["full_dataset_execution_authorized"] is False
    assert [row["cell"] for row in report["attempts"]] == list(CELLS)
    historical = Path(report["inputs"][0]["path"]).parent.parent
    def path(pin):
        result = copied / Path(pin["path"]).relative_to(historical)
        data = result.read_bytes()
        assert len(data) == pin["bytes"] and hashlib.sha256(data).hexdigest() == pin["sha256"]
        return result
    names = {line[1:].split()[0] for pin in report["inputs"] for line in path(pin).read_text().splitlines()
             if line.startswith(">")}
    assert len(names) == 26
    checkpoints, partitions = [], []
    for row in report["attempts"]:
        p, c, r = (row["cell"][i] == "1" for i in (1, 4, 7))
        assert row["status"] == "diagnostic_checked" and row["exit_code"] == 0
        path(row["log"])
        diagnostic = json.loads(path(row["diagnostic"]).read_text())
        assert diagnostic["process_cpu_affinity"] == [0, 1] and diagnostic["cpu"] == 2
        metrics = json.loads(path(diagnostic["metrics"]).read_text())
        expected_stages = {"search", "edge_thresholds", "network_edges", "clustering", "refinement", "orthogroup_materialization"}
        expected_stages |= {"profile_expansion"} if p else set()
        expected_stages |= {"phylogeny_candidates"} if c else set()
        expected_stages |= {"phylogeny"} if r else set()
        assert set(metrics["stages"]) == expected_stages
        assert metrics["metadata"]["native_factorial"] == diagnostic["factors"]
        assert diagnostic["counts"] == row["counts"] == metrics["counts"]
        if p:
            assert row["counts"]["high_sensitivity_profiles"] == 6
            assert row["counts"]["high_sensitivity_profile_hits"] == 26
        if c:
            assert row["counts"]["phylogeny_candidate_merges"] == 0
        if r:
            assert row["counts"]["phylogeny_species_tree_families"] == 5
            assert row["counts"]["phylogeny_reconciled_families"] == 1
            assert row["counts"]["phylogeny_duplications"] == 2
            assert row["counts"]["phylogeny_checkpoint_hits"] == 0
            assert row["counts"]["phylogeny_species_tree_checkpoint_hit"] is False
            phy = copied / row["cell"] / "native/orthohmm_phylogeny"
            assert (phy / "species_tree.rooted.nwk").stat().st_size > 0
            assert (phy / "gene_trees/Family0000005.reconciled.nwk").stat().st_size > 0
            with (phy / "orthohmm_root_hogs.tsv").open() as stream:
                roots = [item["genes"].split(",") for item in csv.DictReader(stream, delimiter="\t")]
            assert len(roots) == 6 and sum(map(len, roots)) == 26
            assert {gene for group in roots for gene in group} == names
        manifest_path = path(row["checkpoint"])
        manifest = json.loads(manifest_path.read_text())
        checkpoints.append(manifest)
        for filename, pin in manifest["files"].items():
            data = (manifest_path.parent / filename).read_bytes()
            assert len(data) == pin["bytes"] and hashlib.sha256(data).hexdigest() == pin["sha256"]
        arrays = [np.load(manifest_path.parent / name, allow_pickle=False) for name in
                  ("gene_to_species.npy", "hit_queries.npy", "hit_targets.npy", "hit_scores.npy")]
        assert len(arrays[0]) == 26 and all(len(array) == 116 for array in arrays[1:])
        assert np.isfinite(arrays[3]).all()
        groups = {tuple(sorted(line.split(":", 1)[1].split()))
                  for line in path(row["native_groups"]).read_text().splitlines()}
        assert len(groups) == 6 and sum(map(len, groups)) == 26
        assert {gene for group in groups for gene in group} == names
        partitions.append(groups)
    assert all(item == checkpoints[0] for item in checkpoints)
    assert all(item == partitions[0] for item in partitions)


def test_all_copied_fixture_artifacts_match_their_inventory():
    root = Path(__file__).resolve().parents[2]
    manifest = json.loads((root / "benchmark_tools/results/native_factorial_adapter_artifacts_20261004.json").read_text())
    copied = root / "benchmark_tools/results/native_factorial_adapter_diagnostic_20261004"
    historical = Path(manifest["roots"][0])
    files = [row for row in manifest["records"] if row["kind"] == "file"]
    assert len(files) == 222 and sum(row["bytes"] for row in files) == 363122
    expected = {Path(row["path"]).relative_to(historical) for row in files}
    assert {path.relative_to(copied) for path in copied.rglob("*") if path.is_file()} == expected
    for row in files:
        data = (copied / Path(row["path"]).relative_to(historical)).read_bytes()
        assert len(data) == row["bytes"] and hashlib.sha256(data).hexdigest() == row["sha256"]

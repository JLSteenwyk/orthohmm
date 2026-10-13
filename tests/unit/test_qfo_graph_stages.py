import csv
import json

import numpy as np
import pytest

from benchmark_tools import trace_qfo_graph_stages as trace


def payload(tmp_path, sources=(0,), targets=(1,), weights=(.8,), dtype=np.int32):
    path = tmp_path / "payload"
    path.mkdir()
    (path / "gene_names.txt").write_text("A\nB\nC\n")
    np.save(path / "sources.npy", np.array(sources, dtype=dtype))
    np.save(path / "targets.npy", np.array(targets, dtype=dtype))
    np.save(path / "weights.npy", np.array(weights, dtype=np.float64))
    return path


def test_native_graph_support_and_weights_are_read_without_reclustering(tmp_path):
    graph = payload(tmp_path, (0, 0, 1), (1, 2, 2), (.8, .7, .9))
    found = trace.watched_graph(graph, ["A", "B", "C"], {(0, 1), (1, 2)},
                                trace.saved_fingerprint(graph), chunk_size=1)
    assert found == {(0, 1): .8, (1, 2): .9}


@pytest.mark.parametrize("sources,targets,weights,match", [
    ((0, 0), (1, 1), (.8, .9), "Duplicate"),
    ((1, 0), (2, 1), (.8, .9), "unordered"),
    ((1,), (0,), (.8,), "endpoints"),
    ((0,), (0,), (.8,), "endpoints"),
    ((-1,), (1,), (.8,), "endpoints"),
    ((0,), (3,), (.8,), "endpoints"),
    ((0,), (1,), (0.,), "weights"),
    ((0,), (1,), (float("nan"),), "Nonfinite")])
def test_invalid_graphs_rejected(tmp_path, sources, targets, weights, match):
    graph = payload(tmp_path, sources, targets, weights)
    with pytest.raises(ValueError, match=match):
        trace.watched_graph(graph, ["A", "B", "C"], {(0, 1)}, trace.saved_fingerprint(graph), chunk_size=1)


def test_graph_dtypes_and_gene_order_are_not_silently_coerced(tmp_path):
    graph = payload(tmp_path, dtype=np.float64)
    with pytest.raises(ValueError, match="dtypes"):
        trace.watched_graph(graph, ["A", "B", "C"], {(0, 1)}, trace.saved_fingerprint(graph))
    with pytest.raises(ValueError, match="gene order"):
        trace.watched_graph(graph, ["B", "A", "C"], {(0, 1)}, trace.saved_fingerprint(graph))


def test_graph_must_reproduce_admitted_native_fingerprint(tmp_path):
    graph = payload(tmp_path)
    wrong = {**trace.saved_fingerprint(graph), "edges": 2}
    with pytest.raises(ValueError, match="fingerprint"):
        trace.watched_graph(graph, ["A", "B", "C"], {(0, 1)}, wrong)


@pytest.mark.parametrize("text,count,match", [
    ("A B\n", 1, "universe"),
    ("A B\nB C\n", 2, "membership"),
    ("A B\nC\n", 3, "group count")])
def test_all_partition_genes_and_recorded_counts_are_checked(tmp_path, text, count, match):
    path = tmp_path / "partition.txt"
    path.write_text(text)
    with pytest.raises(ValueError, match=match):
        trace.partition_memberships(path, {"A", "B", "C"}, {"A", "B"}, count)


def example_rows():
    return [{"left": "alias_A", "right": "alias_B", "reference_label": "f", "category": "FN",
             "stage": "different_candidate_families", "native_gene_left": "A", "native_gene_right": "B",
             "reference_is_cross_species": "True", "same_root_hog": "False", "initial_hit_directions": "2"}]


def stage_data():
    graphs = {name: {(0, 1): .8} for name in trace.GRAPHS}
    partitions = {name: {"A": 0, "B": 0} for name in trace.PARTITIONS}
    for name in ("multipass_refined", "profiles_refined", "candidate", "root"):
        partitions[name] = {"A": 0, "B": 1}
    return graphs, partitions


def test_refinement_transition_is_described_without_changing_pair_truth():
    graphs, partitions = stage_data()
    rows = trace.extend_rows(example_rows(), {"A": 0, "B": 1}, graphs, partitions)
    row = rows[0]
    assert row["category"] == "FN" and row["initial_edge"]
    assert row["same_profile_expanded"] and not row["same_profiles_refined"]
    counts = trace.accounting(rows)
    assert counts["transitions"]["FN/profile_expanded:True/profiles_refined:False"] == 1
    assert counts["stages"]["FN/initial/edge_True"] == 1


@pytest.mark.parametrize("field,value,match", [
    ("same_root_hog", "True", "Root membership"),
    ("stage", "within_candidate_not_predicted", "Candidate membership"),
    ("reference_is_cross_species", "unknown", "boolean"),
    ("category", "FP", "cohort row")])
def test_previous_verified_stage_labels_cannot_be_reinterpreted(field, value, match):
    rows = example_rows()
    rows[0][field] = value
    graphs, partitions = stage_data()
    with pytest.raises(ValueError, match=match):
        trace.extend_rows(rows, {"A": 0, "B": 1}, graphs, partitions)


def test_within_species_rows_are_preserved_in_a_separate_accounting_cohort():
    rows = example_rows()
    rows[0].update(reference_is_cross_species="False", stage="reference_within_species")
    graphs, partitions = stage_data()
    result = trace.extend_rows(rows, {"A": 0, "B": 1}, graphs, partitions)
    assert result[0]["category"] == "FN"
    assert trace.accounting(result)["stages"]["within_species_FN/initial/edge_True"] == 1


def small_run_plan(tmp_path, monkeypatch):
    monkeypatch.setattr(trace, "EXPECTED_GENES", 3)
    graph = payload(tmp_path)
    cohort = tmp_path / "cohort.tsv"
    rows = example_rows()
    with cohort.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    joined, split, root = (tmp_path / name for name in ("joined.txt", "split.txt", "root.tsv"))
    joined.write_text("A B\nC\n")
    split.write_text("A\nB\nC\n")
    root.write_text("root_hog\tsource_family\tgenes\nR0\tF0\tA\nR1\tF0\tB\nR2\tF1\tC\n")
    checked = [trace.record(path) for path in (cohort, joined, split, root, *graph.iterdir())]
    partitions = {name: {"file": trace.record(joined), "groups": 2} for name in trace.PARTITIONS}
    for name in ("multipass_refined", "profiles_refined", "candidate"):
        partitions[name] = {"file": trace.record(split), "groups": 3}
    partitions["root"] = {"file": trace.record(root), "groups": 3}
    plan = {"checked": checked, "sources": [], "source_revision": "fixture", "cohort": trace.record(cohort),
        "checkpoint_names": trace.record(graph / "gene_names.txt"), "expected_reference_pairs": 1,
        "expected_raw_counts": {"TP": 0, "FN": 1}, "partitions": partitions,
        "graphs": {name: {"payload": str(graph), "fingerprint": trace.saved_fingerprint(graph)} for name in trace.GRAPHS}}
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    return path, root


def test_complete_trace_uses_fresh_namespace_and_preserves_inference_outputs(tmp_path, monkeypatch):
    plan, root = small_run_plan(tmp_path, monkeypatch)
    before = root.read_bytes()
    result = trace.run(plan, tmp_path / "result")
    assert result["status"] == "admitted_native_graph_partition_stages_traced"
    assert result["counts"] == {"FN": 1}
    assert root.read_bytes() == before
    assert not result["defaults_changed"] and not result["goal_completed"]
    with pytest.raises(FileExistsError):
        trace.run(plan, tmp_path / "result")


def test_input_mutation_prevents_result_publication(tmp_path, monkeypatch):
    plan, root = small_run_plan(tmp_path, monkeypatch)
    original = trace.extend_rows
    def changed(*args):
        result = original(*args)
        root.write_text("mutated")
        return result
    monkeypatch.setattr(trace, "extend_rows", changed)
    with pytest.raises(ValueError, match="identity changed"):
        trace.run(plan, tmp_path / "result")
    assert not (tmp_path / "result").exists()

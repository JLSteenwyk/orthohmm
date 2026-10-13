import csv
import gzip
import json
import sqlite3

import numpy as np
import pytest

from benchmark_tools import trace_qfo_reference_stages as trace
from orthohmm.accuracy import write_accuracy_checkpoint


def fixture_plan(tmp_path):
    names = [f"sp|{letter}|protein" for letter in "ABCDE"]
    checkpoint = write_accuracy_checkpoint(str(tmp_path / "initial"), names, np.arange(5),
        [0, 1, 1, 3], [1, 0, 2, 4], [1., .9, .8, .7])
    paths = {"checkpoint_manifest": checkpoint / "manifest.json"}
    for key in ("reference", "raw"):
        paths[key] = tmp_path / (key + ".gz")
    with gzip.open(paths["reference"], "wt") as stream:
        stream.write("1\t2\tA\n1\t3\tA\n2\t3\tA\n4\t5\tB\n")
    with gzip.open(paths["raw"], "wt") as stream:
        stream.write("A\tB\tTP\tA\tA\ts0\ts1\n"
                     "A\tC\tFN\tA\tA\ts0\ts2\n"
                     "B\tC\tFN\tA\tA\ts1\ts2\n"
                     "D\tE\tFN\tB\tB\ts3\ts4\n")
    paths["database"] = tmp_path / "participant.db"
    with sqlite3.connect(paths["database"]) as db:
        db.execute("CREATE TABLE proteomes (prot_nr INTEGER, uniprot_id TEXT, species TEXT)")
        db.executemany("INSERT INTO proteomes VALUES (?,?,?)",
                       [(i + 1, letter, f"s{i}") for i, letter in enumerate("ABCDE")])
    paths["candidates"] = tmp_path / "candidates.txt"
    paths["candidates"].write_text(" ".join(names[:3]) + "\n" + "\n".join(names[3:]) + "\n")
    paths["roots"] = tmp_path / "roots.tsv"
    paths["roots"].write_text("root_hog\tsource_family\tgenes\nR0\tF0\t" + ",".join(names[:2]) +
        "\n" + "".join(f"R{i}\tF{i}\t{gene}\n" for i, gene in enumerate(names[2:], start=1)))
    paths["pairs"] = tmp_path / "pairs.tsv"
    paths["pairs"].write_text(f"gene_a\tspecies_a\tgene_b\tspecies_b\n{names[0]}\ts0\t{names[1]}\ts1\n")
    manifest = json.loads(paths["checkpoint_manifest"].read_text())
    for filename in manifest["files"]:
        paths["checkpoint_" + filename] = checkpoint / filename
    arrays = trace.load_accuracy_checkpoint(checkpoint)
    plan = {"inputs": {key: trace.record(path) for key, path in paths.items()},
        "sources": [], "source_revision": "fixture", "expected_universe": trace.check_arrays(*arrays),
        "expected_pairs": 1, "expected_metrics": {"precision": 1., "recall": .25}, "expected_f1": .4}
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    return path, plan, paths


def test_complete_native_trace_preserves_directionality_and_stage_semantics(tmp_path):
    plan, _, _ = fixture_plan(tmp_path)
    result = trace.run(plan, tmp_path / "trace")
    assert result["reference_pairs"] == 4
    assert result["raw_categories"]["counts"] == {"TP": 1, "FP": 0, "FN": 3}
    assert result["summary"] == {
        "TP/predicted/hits_2/root_True": 1,
        "FN/within_candidate_not_predicted/hits_0/root_False": 1,
        "FN/within_candidate_not_predicted/hits_1/root_False": 1,
        "FN/different_candidate_families/hits_1/root_False": 1}
    with (tmp_path / "trace/reference_pair_stages.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert rows[1]["forward_normalized_score"] == "NA"
    assert not result["defaults_promoted"] and not result["goal_completed"]
    with pytest.raises(FileExistsError):
        trace.run(plan, tmp_path / "trace")


def test_root_comembership_is_not_used_as_pair_truth():
    names = ["A", "B"]
    rows = trace.trace_rows({("A", "B"): "family"}, set(), names, {"A": 0, "B": 1},
        np.array([0, 1]), {"A": 0, "B": 0}, {"A": 0, "B": 0}, {})
    assert rows[0]["category"] == "FN" and rows[0]["same_root_hog"]
    assert rows[0]["stage"] == "within_candidate_not_predicted"


def test_reference_row_label_is_preserved_despite_shared_protein_last_label(tmp_path):
    _, _, paths = fixture_plan(tmp_path)
    truth = {(1, 2): "first", (1, 3): "second"}
    mapped = trace.accession_truth(paths["database"], truth, {1: "second", 2: "first", 3: "second"})
    assert mapped == {("A", "B"): "first", ("A", "C"): "second"}


def test_alias_mapping_uses_last_row_without_querying_predictions(tmp_path):
    _, _, paths = fixture_plan(tmp_path)
    with sqlite3.connect(paths["database"]) as db:
        db.execute("INSERT INTO proteomes VALUES (1,'alias','s0')")
    mapped = trace.accession_truth(paths["database"], {(1, 2): "family"}, {1: "family", 2: "family"})
    assert mapped == {("B", "alias"): "family"}


def test_selected_alias_maps_to_exact_database_input_accession(tmp_path):
    plan, manifest, paths = fixture_plan(tmp_path)
    with sqlite3.connect(paths["database"]) as db:
        db.execute("INSERT INTO proteomes VALUES (1,'A_species','s0')")
    with gzip.open(paths["raw"], "rt") as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    with gzip.open(paths["raw"], "wt") as stream:
        for row in rows:
            row[:2] = ["A_species" if value == "A" else value for value in row[:2]]
        csv.writer(stream, delimiter="\t", lineterminator="\n").writerows(rows)
    for key in ("database", "raw"):
        manifest["inputs"][key] = trace.record(paths[key])
    plan.write_text(json.dumps(manifest))
    result = trace.run(plan, tmp_path / "trace")
    assert result["selected_reference_accessions_resolved_via_database_alias"] == 1
    assert result["raw_categories"]["counts"] == {"TP": 1, "FP": 0, "FN": 3}


@pytest.mark.parametrize("index", [{"unrelated": 0}, {"A": 0, "alias": 1}])
def test_missing_or_multiple_input_aliases_are_not_guessed(tmp_path, index):
    _, _, paths = fixture_plan(tmp_path)
    with sqlite3.connect(paths["database"]) as db:
        db.execute("INSERT INTO proteomes VALUES (1,'alias','s0')")
    with pytest.raises(ValueError, match="exactly one input"):
        trace.input_reference_mapping(paths["database"], {1: "f"}, index)


@pytest.mark.parametrize("names", [["sp|A|one", "tr|A|two"], ["A", "A"]])
def test_ambiguous_accession_normalization_rejected(names):
    with pytest.raises(ValueError, match="Ambiguous"):
        trace.normalized_index(names)


def test_missing_or_colliding_database_accessions_rejected(tmp_path):
    _, _, paths = fixture_plan(tmp_path)
    with pytest.raises(ValueError, match="Incomplete"):
        trace.accession_truth(paths["database"], {(1, 6): "family"}, {1: "family", 6: "family"})
    with sqlite3.connect(paths["database"]) as db:
        db.execute("UPDATE proteomes SET uniprot_id='A' WHERE prot_nr=2")
    with pytest.raises(ValueError, match="non-bijective"):
        trace.accession_truth(paths["database"], {(1, 2): "family"}, {1: "family", 2: "family"})


@pytest.mark.parametrize("queries,targets,scores,match", [
    ([0, 0], [1, 1], [1., 2.], "Duplicate"),
    ([0], [5], [1.], "outside"),
    ([0], [1], [float("nan")], "Nonfinite"),
    ([0], [1], [0.], "positive")])
def test_bad_watched_hit_arrays_rejected(queries, targets, scores, match):
    with pytest.raises(ValueError, match=match):
        trace.watched_hits(["A", "B"], np.array([0, 1], dtype=np.int32),
            np.array(queries, dtype=np.int32), np.array(targets, dtype=np.int32),
            np.array(scores, dtype=np.float64), {(0, 1)}, chunk_size=1)


def test_duplicate_unwatched_hits_do_not_claim_global_duplicate_audit():
    hits, _ = trace.watched_hits(["A", "B", "C"], np.arange(3, dtype=np.int32),
        np.array([0, 2, 2], dtype=np.int32), np.array([1, 0, 0], dtype=np.int32),
        np.ones(3), {(0, 1)}, chunk_size=1)
    assert hits == {(0, 1): 1.}


@pytest.mark.parametrize("species,candidates,predicted,match", [
    ([0, 0], {"A": 0, "B": 0}, {("A", "B")}, "within-species"),
    ([0, 1], {"A": 0, "B": 1}, {("A", "B")}, "crosses candidate")])
def test_invalid_reference_or_native_candidate_pair_rejected(species, candidates, predicted, match):
    with pytest.raises(ValueError, match=match):
        trace.trace_rows({("A", "B"): "f"}, predicted, ["A", "B"], {"A": 0, "B": 1},
            np.array(species), candidates, {"A": 0, "B": 1}, {})


def test_within_species_asserted_rows_are_preserved_but_not_recovery_opportunities():
    rows = trace.trace_rows({("A", "B"): "f"}, set(), ["A", "B"], {"A": 0, "B": 1},
        np.array([0, 0]), {"A": 0, "B": 0}, {"A": 0, "B": 0}, {(0, 1): 1.})
    assert len(rows) == 1 and rows[0]["category"] == "FN"
    assert rows[0]["stage"] == "reference_within_species"
    assert not rows[0]["reference_is_cross_species"]
    assert rows[0]["initial_hit_directions"] == 1


def test_native_pair_not_raw_label_controls_trace(tmp_path):
    plan, manifest, paths = fixture_plan(tmp_path)
    names = trace.load_accuracy_checkpoint(paths["checkpoint_manifest"].parent)[0]
    paths["pairs"].write_text(f"gene_a\tspecies_a\tgene_b\tspecies_b\n{names[0]}\ts0\t{names[2]}\ts2\n")
    manifest["inputs"]["pairs"] = trace.record(paths["pairs"])
    plan.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Native truth pairs differ"):
        trace.run(plan, tmp_path / "trace")
    assert not (tmp_path / "trace").exists()


@pytest.mark.parametrize("kind", ["candidates", "roots"])
def test_partition_gene_universe_must_be_complete(tmp_path, kind):
    plan, manifest, paths = fixture_plan(tmp_path)
    paths[kind].write_text("\n".join(paths[kind].read_text().splitlines()[:-1]) + "\n")
    manifest["inputs"][kind] = trace.record(paths[kind])
    plan.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="universe"):
        trace.run(plan, tmp_path / "trace")


def test_evidence_mutation_prevents_result_publication(tmp_path, monkeypatch):
    plan, _, paths = fixture_plan(tmp_path)
    original = trace.trace_rows
    def changed(*args):
        rows = original(*args)
        paths["roots"].write_text("mutated")
        return rows
    monkeypatch.setattr(trace, "trace_rows", changed)
    with pytest.raises(ValueError, match="identity changed"):
        trace.run(plan, tmp_path / "trace")
    assert not (tmp_path / "trace").exists()


@pytest.mark.parametrize("field,value,match", [
    ("expected_pairs", 2, "pair count"),
    ("expected_f1", .5, "metric differs"),
    ("expected_universe", {}, "universe differs")])
def test_admitted_counts_and_scores_cannot_be_silently_changed(tmp_path, field, value, match):
    plan, manifest, _ = fixture_plan(tmp_path)
    manifest[field] = value
    plan.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match=match):
        trace.run(plan, tmp_path / "trace")
    assert not (tmp_path / "trace").exists()


def test_duplicate_native_pairs_rejected(tmp_path):
    _, _, paths = fixture_plan(tmp_path)
    with paths["pairs"].open("a") as stream:
        stream.write(paths["pairs"].read_text().splitlines()[1] + "\n")
    with pytest.raises(ValueError, match="duplicate"):
        trace.native_truth_pairs(paths["pairs"], trace.load_accuracy_checkpoint(paths["checkpoint_manifest"].parent)[0],
                                {("A", "B")}, 2)


def test_checksum_mismatch_rejected_before_any_results(tmp_path):
    plan, _, paths = fixture_plan(tmp_path)
    paths["pairs"].write_text("altered")
    with pytest.raises(ValueError, match="identity changed"):
        trace.run(plan, tmp_path / "trace")
    assert not (tmp_path / "trace").exists()

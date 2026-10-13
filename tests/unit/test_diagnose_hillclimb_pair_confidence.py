import gzip
import json

import pytest

from benchmark_tools.diagnose_hillclimb_pair_confidence import HEADER, read_confidence, run, summarize
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_membership_hillclimb import PRODUCTS


def native(tmp_path, rows):
    path = tmp_path / "confidence.tsv"
    path.write_text("\t".join(HEADER) + "\n" + "".join("\t".join(row) + "\n" for row in rows))
    return path


def test_normalized_scored_pairs_have_exact_native_confidence_coverage(tmp_path):
    path = native(tmp_path, [("sp|A|A", "s1", "tr|B|B", "s2", "high"),
                            ("sp|A|A", "s1", "tr|C|C", "s3", "medium")])
    found, counts = read_confidence(path, {("A", "B")}, 2)
    assert found == {("A", "B"): "high"} and counts == {"high": 1, "medium": 1}


@pytest.mark.parametrize("mutation", ["header", "same_species", "same_gene", "unknown", "duplicate", "missing", "count", "alias"])
def test_bad_confidence_tables_are_rejected(tmp_path, mutation):
    rows = [("sp|A|A", "s1", "tr|B|B", "s2", "high")]
    watched, count = {("A", "B")}, 1
    if mutation == "same_species":
        rows[0] = ("sp|A|A", "s1", "tr|B|B", "s1", "high")
    elif mutation == "same_gene":
        rows[0] = ("sp|A|A", "s1", "sp|A|A", "s2", "high")
    elif mutation == "unknown":
        rows[0] = (*rows[0][:4], "low")
    elif mutation == "duplicate":
        rows *= 2
        count = 2
    elif mutation == "missing":
        watched.add(("C", "D"))
    elif mutation == "count":
        count = 2
    elif mutation == "alias":
        rows.append(("sp|A|OTHER", "s1", "tr|C|C", "s2", "high"))
        count = 2
    path = native(tmp_path, rows)
    if mutation == "header":
        path.write_text("wrong\n")
    with pytest.raises(ValueError):
        read_confidence(path, watched, count)


def test_categories_keep_removed_added_and_shared_confidence_transitions_separate():
    left = {("A", "B"): 500000, ("C", "D"): 900000}
    right = {("A", "B"): 500000, ("E", "F"): 300000}
    result = summarize(left, right, {("A", "B"): "high", ("C", "D"): "medium"},
                       {("A", "B"): "medium", ("E", "F"): "high"})
    assert [(r["category"], r["confidence"], r["scored_relations"], r["serialized_mean"]) for r in result] == [
        ("added", "high", 1, .3), ("removed", "medium", 1, .9), ("shared", "high->medium", 1, .5)]


def test_changed_shared_annotation_score_is_not_a_pair_set_only_diagnostic():
    with pytest.raises(ValueError, match="Shared serialized"):
        summarize({("A", "B"): 500000}, {("A", "B"): 400000}, {("A", "B"): "high"}, {("A", "B"): "high"})


def fixture_run(tmp_path):
    cases, cells = [], []
    for label in ("control", "deferred500"):
        root = tmp_path / label
        root.mkdir()
        confidence = native(root, [("sp|A|A", "s1", "tr|B|B", "s2", "high")])
        for metric in ("GO", "EC"):
            directory = root / "qfo_scoring/results" / metric
            directory.mkdir(parents=True)
            with gzip.open(directory / "endpoint_raw.txt.gz", "wt") as stream:
                stream.write(f"# {metric} Similarities between orthologs from fixture\n# Computing timestamp: test\n"
                             f"# Protein ID 1<tab>Protein ID 2<tab>{metric} Similarity\nA\tB\t0.500000\n")
        cases.append({"label": label, "output": str(root)})
        cells.append({"method_key": "deferred_inferred:" + label, "dataset": "QfO", "status": "independent_content_readback_verified",
            "scores": {"QfO: GO": .5, "QfO: EC": .5}, "details": {"GO": {"assessed_relations": 1}, "EC": {"assessed_relations": 1}},
            "products": {PRODUCTS[2]: record(confidence)}, "coverage": {"pairs": 1}})
    plan = tmp_path / "execution_plan.json"
    plan.write_text(json.dumps({"cases": cases}))
    full = tmp_path / "full.json"
    full.write_text(json.dumps({"cells": cells, "captures": [{"source_path": str(plan), "snapshot": record(plan)}]}))
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps({"full_report": record(full)}))
    return receipt, full, cells


def test_receipt_bound_diagnostic_preserves_original_files_and_scores(tmp_path):
    receipt, _, _ = fixture_run(tmp_path)
    before = {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
    result = run(receipt, tmp_path / "diagnostic.json")
    assert set(result["comparisons"]) == {"GO", "EC"}
    assert result["comparisons"]["GO"][0]["confidence"] == "high->high"
    assert all(p.read_bytes() == data for p, data in before.items())
    with pytest.raises(ValueError, match="existing"):
        run(receipt, tmp_path / "diagnostic.json")


@pytest.mark.parametrize("mutation", ["checksum", "order", "status", "mean", "count"])
def test_unadmitted_or_inconsistent_raw_endpoints_fail_closed(tmp_path, mutation):
    receipt, full, cells = fixture_run(tmp_path)
    payload = json.loads(full.read_text())
    if mutation == "checksum":
        payload["tampered"] = True
    elif mutation == "order":
        payload["cells"].reverse()
    elif mutation == "status":
        payload["cells"][0]["status"] = "native_running"
    elif mutation == "mean":
        payload["cells"][0]["scores"]["QfO: GO"] = .9
    elif mutation == "count":
        payload["cells"][0]["details"]["GO"]["assessed_relations"] = 2
    full.write_text(json.dumps(payload))
    if mutation != "checksum":
        receipt.write_text(json.dumps({"full_report": record(full)}))
    with pytest.raises(ValueError):
        run(receipt, tmp_path / "diagnostic.json")

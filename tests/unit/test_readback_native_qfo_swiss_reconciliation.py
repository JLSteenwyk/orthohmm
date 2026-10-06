"""Actual saved-tree readback and corruption/refusal checks."""

import copy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import readback_native_qfo_swiss_reconciliation as reader

RESULTS = Path(reader.__file__).parent / "results"
REPORT = RESULTS / "native_qfo_swiss_reconciliation_trace_20261006_v1.json"
REPORT_SHA = "53c692cf7dab9affac7b0240cf5191d2530bd4369b66f4165e8acb811391dacb"


def test_actual_independent_newick_readback():
    result = reader.verify(REPORT, REPORT_SHA)
    assert (result["source_families_checked"], result["tree_leaves_checked"],
            result["distinct_exclusion_lcas_checked"], result["changed_pairs_checked"]) == (23, 1139, 51, 2023)
    assert result == json.loads((RESULTS / "native_qfo_swiss_reconciliation_newick_readback_20261006.json").read_text())
    assert result["original_node_annotation_admission_established"] is False


@pytest.mark.parametrize("fault", ("source", "scope", "admission_scope", "summary", "inventory", "lca", "endpoint",
                                 "root_marker", "pair_event", "label", "missing_pair"))
def test_corrupt_localization_refuses(tmp_path, fault):
    report = copy.deepcopy(json.loads(REPORT.read_text()))
    if fault == "source":
        report["source"] = reader.record(__file__)
    elif fault == "scope":
        report["uncertainty_admitted"] = True
    elif fault == "admission_scope":
        report["node_annotations_previously_inventoried"] = True
    elif fault == "summary":
        report["summary"]["TP_same_root"] += 1
    elif fault == "inventory":
        report["selected_source_families"] += 1
    else:
        with open(report["pair_ledger"]["path"], newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        if fault == "lca":
            rows[0]["lca_node"] = "G99999"
        elif fault == "endpoint":
            rows[0]["gene_a"] = "missing"
        elif fault == "root_marker":
            rows[0]["same_root_hog"] = "False"
        elif fault == "pair_event":
            rows[0]["pair_event"] = "uncertain"
        elif fault == "label":
            rows[0]["before"] = "TP"
        elif fault == "missing_pair":
            rows.pop()
        path = tmp_path / "pairs.tsv"
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        report["pair_ledger"] = reader.record(path)
    path = tmp_path / "report.json"
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError):
        reader.verify(path, reader.record(path)["sha256"])


def test_changed_report_digest_refuses():
    with pytest.raises(ValueError, match="Changed localization report"):
        reader.verify(REPORT, "0" * 64)


@pytest.mark.parametrize("fault", (None, "family", "species", "status", "hash", "leaves", "topology"))
def test_checkpoint_tree_identity(tmp_path, fault):
    directory = tmp_path
    (directory / "checkpoints").mkdir()
    (directory / "gene_trees").mkdir()
    family = "Family0000000"
    raw, rooted = "((a,b),(c,d));\n", "[&R] ((a,b),(c,d));\n"
    annotated = "[&R] ((a,b)'G00000|S@S0',(c,d)'G00001|S@S0')'G00002|D@S0';\n"
    if fault == "leaves":
        annotated = annotated.replace("a,b", "a,x")
    elif fault == "topology":
        annotated = annotated.replace("a,b", "a,c").replace("c,d", "b,d")
    checkpoint = dict(schema_version=2, status="complete", family_id=family,
                      species_tree_sha256="species", genes=list("abcd"))
    for suffix, key, text in (("raw", "raw_tree_sha256", raw), ("rooted", "rooted_tree_sha256", rooted),
                             ("reconciled", "annotated_tree_sha256", annotated)):
        path = directory / "gene_trees" / (family + "." + suffix + ".nwk")
        path.write_text(text)
        checkpoint[key] = reader.record(path)["sha256"]
    if fault == "family":
        checkpoint["family_id"] = "Family0000001"
    elif fault == "species":
        checkpoint["species_tree_sha256"] = "other"
    elif fault == "status":
        checkpoint["status"] = "running"
    elif fault == "hash":
        checkpoint["raw_tree_sha256"] = "0" * 64
    (directory / "checkpoints" / (family + ".json")).write_text(json.dumps(checkpoint))
    if fault is None:
        evidence = []
        tree = reader.checked_tree(directory, family, "species", evidence)
        assert len(tree.get_terminals()) == len(evidence) == 4
    else:
        with pytest.raises(ValueError):
            reader.checked_tree(directory, family, "species", [])


@pytest.mark.parametrize("family", ("../escape", "Family1", "Family0000000/else"))
def test_invalid_family_path_refuses(tmp_path, family):
    with pytest.raises(ValueError, match="family ID"):
        reader.checked_tree(tmp_path, family, "species", [])

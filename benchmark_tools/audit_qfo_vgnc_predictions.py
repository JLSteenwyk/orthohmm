"""Reconstruct VGNC categories from prediction databases and compare exact pairs."""

import argparse
from collections import defaultdict
import gzip
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_vgnc_mapping import reference_data
from benchmark_tools.report_qfo_recovered_stages import ADMISSION_SHA, STAGES
from benchmark_tools.snapshot_orthohmm_input_order import record


def classify(truth, predictions, labels, species):
    families = defaultdict(set)
    for protein, family in labels.items():
        families[species[protein]].add(family)
    # FP eligibility is independent of asserted truth in the native scorer.
    # Do not silently force disjoint TP/FP sets when shared labels permit overlap.
    fp = {(a, b) for a, b in predictions if labels[a] != labels[b]
          and labels[a] in families[species[b]] and labels[b] in families[species[a]]}
    truth = set(truth)
    return {"TP": truth & predictions, "FN": truth - predictions, "FP": fp}


def pair_digest(pairs):
    return hashlib.sha256(json.dumps(sorted(pairs), separators=(",", ":")).encode()).hexdigest()


def read_predictions(database, labels):
    metadata, predictions = {}, set()
    with sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True) as db:
        for protein, accession, species in db.execute(
            "SELECT prot_nr, uniprot_id, species FROM proteomes ORDER BY rowid"
        ):
            if protein in labels:
                if protein in metadata and metadata[protein][1] != species:
                    raise ValueError("Conflicting species mapping")
                metadata[protein] = (accession, species)
        if set(metadata) != set(labels) or len({a for a, _ in metadata.values()}) != len(labels):
            raise ValueError("Incomplete or non-bijective final accession mapping")
        for protein in sorted(labels):
            for a, b in db.execute(
                "SELECT prot_nr1, prot_nr2 FROM orthologs WHERE prot_nr1 = ? AND prot_nr2 > ?",
                (protein, protein),
            ):
                if b in labels:
                    predictions.add((a, b))
    return metadata, predictions


def compare_raw(raw, expected, metadata):
    by_accession = {accession: protein for protein, (accession, _) in metadata.items()}
    observed = {label: set() for label in ("TP", "FP", "FN")}
    with gzip.open(raw, "rt") as stream:
        for line in stream:
            a, b, label, *_ = line.rstrip("\n").split("\t")
            pair = tuple(sorted((by_accession[a], by_accession[b])))
            if label not in observed or pair in observed[label]:
                raise ValueError("Invalid or duplicate raw category")
            observed[label].add(pair)
    for label in observed:
        if observed[label] != expected[label]:
            missing = len(expected[label] - observed[label])
            extra = len(observed[label] - expected[label])
            raise ValueError(f"Raw {label} mismatch: {missing} missing, {extra} extra")


def audit(repo):
    results = repo / "benchmark_tools/results"
    mapping_path = results / "vgnc_mapping_audit_20260917.json"
    mapping = json.loads(mapping_path.read_text())
    inventory_path = results / "vgnc_raw_label_inventory_20260917.json"
    inventory = json.loads(inventory_path.read_text())
    reference_inventory = json.loads(Path(inventory["reference_inventory"]["path"]).read_text())
    reference = reference_inventory["reference"]
    admission_path = results / "qfo_recovered_assessment_admitted_20260917.json"
    admission_identity = record(admission_path)
    if admission_identity["sha256"] != ADMISSION_SHA:
        raise ValueError("Changed admitted endpoints")
    admission = json.loads(admission_path.read_text())
    if [r["stage"] for r in admission["records"]] != list(STAGES):
        raise ValueError("Changed stage ordering")
    checked = [record(mapping_path), record(inventory_path), inventory["reference_inventory"],
               reference, admission_identity, inventory["native_scorer"]]
    truth, labels = reference_data(Path(reference["path"]))
    stages = []
    for index, item in enumerate(inventory["records"]):
        if item["stage_index"] != index or mapping["stages"][index]["stage_index"] != index:
            raise ValueError("Changed inventory stage ordering")
        database = Path(mapping["stages"][index]["database"])
        identity = record(database)
        checked.extend((identity, item["raw"]))
        metadata, predictions = read_predictions(database, labels)
        categories = classify(truth, predictions, labels, {p: s for p, (_, s) in metadata.items()})
        compare_raw(Path(item["raw"]["path"]), categories, metadata)
        counts = {label: len(pairs) for label, pairs in categories.items()}
        tp, fp, fn = (counts[label] for label in ("TP", "FP", "FN"))
        scores = {"PPV": tp / (tp + fp), "TPR": tp / (tp + fn), "F1": 2 * tp / (2 * tp + fp + fn)}
        endpoint = admission["records"][index]["assessment"]["endpoints"]["VGNC"]
        for metric, axis in (("TPR", "metric_x"), ("PPV", "metric_y")):
            if not math.isclose(scores[metric], endpoint["native_participant"][axis], rel_tol=0, abs_tol=5e-8):
                raise ValueError("Native metric mismatch: " + metric)
        if not math.isclose(scores["F1"], endpoint["score"], rel_tol=0, abs_tol=5e-8):
            raise ValueError("Native harmonic F1 mismatch")
        stages.append({"stage": STAGES[index], "database": identity, "counts": counts, "scores": scores,
                       "predictions_among_reference_proteins": len(predictions),
                       "unscored_predictions": len(predictions - categories["TP"] - categories["FP"]),
                       "tp_fp_overlap": len(categories["TP"] & categories["FP"]),
                       "prediction_pairs_sha256": pair_digest(predictions),
                       "category_pairs_sha256": {label: pair_digest(pairs) for label, pairs in categories.items()}})
    for identity in checked:
        if record(identity["path"]) != identity:
            raise ValueError("Evidence changed during audit: " + identity["path"])
    return {"status": "exact_vgnc_prediction_categories_and_scores_verified",
            "source": record(__file__), "reference_reader": record(Path(__file__).with_name("audit_qfo_vgnc_mapping.py")),
            "checked_inputs": checked, "stages": stages,
            "limitations": ["Four recovered stages only, not all publication competitors.",
                            "Validates scoring of mapped prediction databases, not upstream pair conversion completeness.",
                            "No independent-family uncertainty or new biological truth established.",
                            "Retains native one-direction query and last-row alias behavior."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(audit(args.repo), indent=2, sort_keys=True) + "\n")

"""Audit VGNC reference/accession mapping without assuming independent families."""

import argparse
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.snapshot_orthohmm_input_order import record


def reference_data(path):
    truth, labels = {}, {}
    with gzip.open(path, "rt") as stream:
        for line in stream:
            left, right, family = line.rstrip("\n").split("\t")
            pair = tuple(sorted((int(left), int(right))))
            if pair[0] == pair[1] or pair in truth:
                raise ValueError("Duplicate or self reference pair")
            truth[pair] = family
            # Preserve the native scorer's ordered last-label assignment.
            for protein in pair:
                labels[protein] = family
    return truth, labels


def mapped_reference(database, truth, labels):
    mapped, accession_owners = {}, {}
    extra_rows = {"identical": 0, "alias": 0}
    with sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True) as connection:
        for protein, accession, species in connection.execute(
            "SELECT prot_nr, uniprot_id, species FROM proteomes ORDER BY rowid"
        ):
            if protein not in labels:
                continue
            if mapped.get(protein) == (accession, labels[protein], species):
                extra_rows["identical"] += 1
                continue
            if accession in accession_owners and accession_owners[accession] != protein:
                raise ValueError("Non-bijective reference accession mapping")
            if protein in mapped:
                if mapped[protein][2] != species:
                    raise ValueError("Conflicting species for reference protein")
                extra_rows["alias"] += 1
            # Native dictionary assignment keeps the last encountered alias.
            # Raw comparison below tests this reconstruction on retained outputs.
            mapped[protein] = (accession, labels[protein], species)
            accession_owners[accession] = protein
    if set(mapped) != set(labels):
        raise ValueError("Missing reference proteins in database")
    annotations = {value[0]: value[1:] for value in mapped.values()}
    pairs = {tuple(sorted((mapped[a][0], mapped[b][0]))) for a, b in truth}
    digest = hashlib.sha256(json.dumps(sorted(mapped.items()), separators=(",", ":")).encode()).hexdigest()
    return pairs, annotations, digest, extra_rows


def validate_raw(path, truth, annotations):
    observed = {label: set() for label in ("TP", "FP", "FN")}
    species_families = defaultdict(set)
    for family, species in annotations.values():
        species_families[species].add(family)
    with gzip.open(path, "rt") as stream:
        for line in stream:
            a, b, label, fa, fb, sa, sb = line.rstrip("\n").split("\t")
            if annotations.get(a) != (fa, sa) or annotations.get(b) != (fb, sb):
                raise ValueError("Raw protein annotation differs from reference/database")
            pair = tuple(sorted((a, b)))
            if label not in observed or a == b or pair in observed[label]:
                raise ValueError("Invalid category, self pair or duplicate raw row")
            if label in ("TP", "FN") and pair not in truth:
                raise ValueError("Raw truth label not asserted by reference")
            if label == "FP" and not (
                fa != fb and fa in species_families[sb] and fb in species_families[sa]
            ):
                raise ValueError("Raw false positive fails native eligibility rule")
            observed[label].add(pair)
    if observed["TP"] & observed["FN"] or observed["TP"] | observed["FN"] != truth:
        raise ValueError("Raw truth partition differs from mapped reference")
    return {"counts": {key: len(value) for key, value in observed.items()},
            "tp_fp_overlap": len(observed["TP"] & observed["FP"]),
            "fn_fp_overlap": len(observed["FN"] & observed["FP"])}


def audit(repo):
    results = repo / "benchmark_tools/results"
    inventory_path = results / "vgnc_raw_label_inventory_20260917.json"
    inventory = json.loads(inventory_path.read_text())
    ref_inventory = json.loads(Path(inventory["reference_inventory"]["path"]).read_text())
    checked = [record(inventory_path), inventory["reference_inventory"],
               ref_inventory["reference"], inventory["native_scorer"]]
    truth, labels = reference_data(Path(ref_inventory["reference"]["path"]))
    stages, common_mapping = [], None
    for item in inventory["records"]:
        index = item["stage_index"]
        database = repo / f"qfo_benchmark/scoring/checked_v2_{index}/other/ohmm_checked_v2_{index}.db"
        before = database.stat()
        pairs, annotations, digest, duplicates = mapped_reference(database, truth, labels)
        if common_mapping is not None and digest != common_mapping:
            raise ValueError("Reference mapping differs across stages")
        common_mapping = digest
        checked.append(item["raw"])
        validation = validate_raw(Path(item["raw"]["path"]), pairs, annotations)
        if validation["counts"] != item["counts"]:
            raise ValueError("Raw counts differ from retained inventory")
        after = database.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError("Database changed during mapping audit")
        stages.append({"stage_index": index, "database": str(database),
                       "database_bytes": before.st_size,
                       "extra_reference_mapping_rows": duplicates,
                       "selected_reference_mapping_sha256": digest, **validation})
    for item in checked:
        if record(item["path"]) != item:
            raise ValueError("Evidence identity changed: " + item["path"])
    return {"status": "reference_mapping_and_raw_labels_verified", "source": record(__file__),
            "checked_inputs": checked, "reference_pairs": len(truth),
            "reference_proteins": len(labels), "stages": stages,
            "limitations": ["Selected mapping digest is not a checksum of the full prediction database.",
                            "Last rowid alias reconstructs retained raw labels; native SQL has no explicit ORDER BY, so no general alias-order guarantee is implied.",
                            "Does not independently rescore database predictions or verify omitted false positives.",
                            "Does not establish independent family units or confidence intervals."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(audit(args.repo), indent=2, sort_keys=True) + "\n")

"""Compare archive sequence content to the retained scorer, without extraction."""

import argparse
import gzip
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.compare_qfo_corrected_archive import compare, missing_by_species
from benchmark_tools.audit_qfo_input_sequences import compare_database
from benchmark_tools.inventory_swiss_sequences import PREPARED_SHA
from benchmark_tools.resolve_swiss_sequence_aliases import MAPPING_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

DATABASE_SHA = "f71e282f504d007306a34da98ab717c2c7fad0a28e4d503caa47a042a41ec9a8"


def classify_matches(result):
    explained, unexplained = [], []
    for row in result["differences"]:
        (explained if row["BOUZ_to_X_sequence"] == row["native_sequence"] else unexplained).append(row["numeric_id"])
    if result["sequence_identical"] + len(explained) + len(unexplained) != result["mapped_input_sequences"]:
        raise ValueError("Sequence comparison counts inconsistent")
    return {"exact_sequence_matches": result["sequence_identical"],
            "BOUZ_to_X_only_numeric_ids": explained, "unexplained_sequence_difference_ids": unexplained,
            "all_mapped_sequences_accounted_for_by_exact_or_BOUZ_to_X": not unexplained}


def audit(archive, prepared_path, mapping_path, database):
    sources = [record(p) for p in (archive, prepared_path, mapping_path, database)]
    if [r["sha256"] for r in sources[1:]] != [PREPARED_SHA, MAPPING_SHA, DATABASE_SHA]:
        raise ValueError("Changed frozen reference source")
    prepared = json.loads(prepared_path.read_text())
    with gzip.open(mapping_path, "rt") as stream:
        mapping = json.load(stream)
    expected = {Path(r["path"]).name: r for r in prepared["input_fastas"]}
    if len(expected) != 78:
        raise ValueError("Expected 78 original proteomes")
    sequences = {}
    archive_result = compare(archive, expected, mapping["mapping"], set(), sequence_records=sequences)
    archive_result["missing_numeric_ids_by_species"] = missing_by_species(
        archive_result["mapping_numeric_ids_without_canonical_accession"], mapping["species"], mapping["Goff"])
    native = compare_database(database, sequences, mapping["Goff"][-1])
    classes = classify_matches(native)
    for identity in sources:
        check(identity)
    return {"status": "archive_native_sequence_comparison", "inputs": sources, "source": record(__file__),
            "helpers": [record(Path(__file__).with_name(name)) for name in
                        ("compare_qfo_corrected_archive.py", "audit_qfo_input_sequences.py")],
            "archive_comparison": archive_result, "native_sequence_comparison": native,
            "sequence_match_classes": classes,
            "complete_mapping_coverage": not archive_result["mapping_numeric_ids_without_canonical_accession"] and archive_result["unmapped_accessions"] == 0,
            "inputs_modified": False,
            "limitations": ["BOUZ-to-X is explicitly reported representation equivalence, not exact input-byte equality.",
                            "No input extracted or normalized on disk; no inference, accuracy evaluation or default change performed.",
                            "Sequence/mapping compatibility does not independently validate annotations, publisher identity or every scorer implementation."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("archive", "prepared", "mapping", "database", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.archive, args.prepared, args.mapping, args.database)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")

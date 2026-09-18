"""Describe exact residue differences without normalizing benchmark inputs."""

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_input_sequences import parse_entry, sequence_identity
from benchmark_tools.inventory_swiss_sequences import AMINO_ACIDS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

AUDIT_SHA = "480a09d0c60c95274d93fb49e15ae95b45a99d0cdcd17c72bf9b3af066cc12c4"


def classify(original, native):
    if original == native:
        raise ValueError("Expected a verified difference")
    if len(original) != len(native):
        return {"category": "different_length_unaligned", "length_delta_native_minus_input": len(native) - len(original),
                "residue_changes": None, "changed_positions": None}
    changes = Counter((a, b) for a, b in zip(original, native) if a != b)
    only_noncanonical_to_x = all(a not in AMINO_ACIDS and b == "X" for a, b in changes)
    category = "noncanonical_to_X_only" if only_noncanonical_to_x else "other_same_length"
    return {"category": category, "length_delta_native_minus_input": 0,
            "residue_changes": [{"input": a, "native": b, "count": n} for (a, b), n in sorted(changes.items())],
            "changed_positions": sum(changes.values())}


def audit(audit_path):
    identity = record(audit_path)
    if identity["sha256"] != AUDIT_SHA:
        raise ValueError("Changed original sequence audit")
    previous = json.loads(audit_path.read_text())
    rows = previous["differences"]
    wanted = {r["accession"]: r for r in rows}
    by_number = {r["numeric_id"]: r for r in rows}
    if len(wanted) != len(rows) or len(by_number) != len(rows):
        raise ValueError("Duplicate difference identity")
    filenames = {r["source_file"] for r in rows}
    inputs = [r for r in previous["fasta_inputs"] if Path(r["path"]).name in filenames]
    if len(inputs) != len(filenames):
        raise ValueError("Missing source FASTAs")
    sequences = {}
    for item in inputs:
        check(item)
        for entry in SeqIO.parse(item["path"], "fasta"):
            accession = entry.id.split("|")[1]
            if accession not in wanted:
                continue
            expected = wanted[accession]
            if accession in sequences or Path(item["path"]).name != expected["source_file"]:
                raise ValueError("Ambiguous selected sequence")
            sequence = str(entry.seq)
            if sequence_identity(sequence) != expected["sequence"]:
                raise ValueError("Input sequence identity changed")
            sequences[accession] = sequence
    if set(sequences) != set(wanted):
        raise ValueError("Missing selected input sequence")
    database = previous["inputs"][2]
    check(database)
    results = []
    with Path(database["path"]).open() as stream:
        for number, line in enumerate(stream, 1):
            if number not in by_number:
                continue
            expected = by_number[number]
            species, aliases, native_identity = parse_entry(line)
            if species != expected["species"] or expected["accession"] not in aliases or native_identity != expected["native_sequence"]:
                raise ValueError("Native sequence provenance changed")
            native = ET.fromstring(line).findtext("SEQ")
            results.append({**expected, **classify(sequences[expected["accession"]], native)})
    if len(results) != len(rows):
        raise ValueError("Missing native difference entries")
    categories, per_species, substitutions = Counter(), defaultdict(Counter), Counter()
    for row in results:
        categories[row["category"]] += 1
        per_species[row["species"]][row["category"]] += 1
        for change in row["residue_changes"] or []:
            substitutions[(change["input"], change["native"])] += change["count"]
    for item in [identity, database, *inputs]:
        check(item)
    return {"status": "exact_sequence_differences_classified", "source": record(__file__),
            "reader": record(Path(__file__).with_name("audit_qfo_input_sequences.py")),
            "input_audit": identity, "database": database, "selected_fastas": inputs,
            "categories": dict(categories), "species_categories": {k: dict(v) for k, v in per_species.items()},
            "same_length_substitutions": [{"input": a, "native": b, "positions": n} for (a, b), n in sorted(substitutions.items())],
            "differences": results, "inputs_normalized": False,
            "limitations": ["Exact positional comparison only for equal-length sequences; different lengths are not aligned.",
                            "Observed symbol replacement does not itself establish the historical preprocessing implementation.",
                            "No benchmark inputs, reference identities or scores modified; corrected-release sequence verification remains separate."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.audit)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")

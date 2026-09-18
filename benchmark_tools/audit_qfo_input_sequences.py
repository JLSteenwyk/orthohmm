"""Compare frozen input sequence bytes with mapped native scorer database entries."""

import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.inventory_swiss_sequences import PREPARED_SHA
from benchmark_tools.resolve_swiss_sequence_aliases import MAPPING_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def sequence_identity(sequence):
    if not sequence or any(c.isspace() for c in sequence):
        raise ValueError("Empty or whitespace-containing sequence")
    return {"length": len(sequence), "sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest()}


def parse_entry(line):
    entry = ET.fromstring(line)
    if entry.tag != "E":
        raise ValueError("Unexpected native entry root")
    values = {}
    for tag in ("OS", "MAPIDS", "SEQ"):
        nodes = entry.findall(tag)
        if len(nodes) != 1 or len(nodes[0]) or not nodes[0].text:
            raise ValueError("Missing, duplicate or nested native field")
        values[tag] = nodes[0].text
    aliases = values["MAPIDS"].split("; ")
    if any(not value for value in aliases) or len(set(aliases)) != len(aliases):
        raise ValueError("Empty or duplicate native alias")
    return values["OS"], set(aliases), sequence_identity(values["SEQ"])


def compare_database(path, mapped_inputs, expected_entries):
    counts, mismatches, missing = Counter(), [], Counter()
    total = 0
    with path.open() as stream:
        for number, line in enumerate(stream, 1):
            species, aliases, sequence = parse_entry(line)
            total = number
            if number not in mapped_inputs:
                missing[species] += 1
                continue
            item = mapped_inputs[number]
            if item["accession"] not in aliases:
                raise ValueError("Input mapping disagrees with native entry alias")
            same = sequence == item["sequence"]
            counts["sequence_identical" if same else "sequence_different"] += 1
            if not same:
                mismatches.append({"numeric_id": number, "species": species, **item, "native_sequence": sequence})
    if total != expected_entries or counts.total() != len(mapped_inputs):
        raise ValueError("Native entry count or mapped input coverage differs")
    return {"native_entries": total, "mapped_input_sequences": len(mapped_inputs),
            "sequence_identical": counts["sequence_identical"], "sequence_different": counts["sequence_different"],
            "differences": mismatches, "native_entries_without_input_by_species": dict(missing)}


def audit(prepared_path, mapping_path, database, database_sha):
    sources = [record(p) for p in (prepared_path, mapping_path, database)]
    if [r["sha256"] for r in sources] != [PREPARED_SHA, MAPPING_SHA, database_sha]:
        raise ValueError("Changed frozen sequence or mapping source")
    prepared = json.loads(prepared_path.read_text())
    with gzip.open(mapping_path, "rt") as stream:
        mapping = json.load(stream)
    matched, accessions, unmapped = {}, set(), []
    for identity in prepared["input_fastas"]:
        check(identity)
        for entry in SeqIO.parse(identity["path"], "fasta"):
            fields = entry.id.split("|")
            if len(fields) != 3 or fields[0] not in {"sp", "tr"} or not all(fields):
                raise ValueError("Unexpected input header")
            accession = fields[1]
            if accession in accessions:
                raise ValueError("Duplicate input accession")
            accessions.add(accession)
            number = mapping["mapping"].get(accession)
            if number is None:
                unmapped.append(accession)
                continue
            if type(number) is not int or number <= 0 or number in matched:
                raise ValueError("Invalid or nonunique input numeric identity")
            matched[number] = {"accession": accession, "source_file": Path(identity["path"]).name,
                               "sequence": sequence_identity(str(entry.seq))}
        check(identity)
    result = compare_database(database, matched, mapping["Goff"][-1])
    for identity in sources:
        check(identity)
    return {"status": "original_input_native_sequence_comparison", **result,
            "input_accessions": len(accessions), "unmapped_accessions": sorted(unmapped),
            "inputs": sources, "fasta_inputs": prepared["input_fastas"], "source": record(__file__),
            "limitations": ["Exact case-sensitive sequence bytes; no alignment or normalization masks substitutions.",
                            "Native database entry ordinals are checked against mapped input aliases and total reference count.",
                            "XML record parsing is independent of Darwin; does not validate annotation correctness or every scorer algorithm.",
                            "This compares original frozen inputs only; corrected release and historical input parity remain separate."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("prepared", "mapping", "database", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--database-sha256", required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.prepared, args.mapping, args.database, args.database_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")

"""Compare retained OrthoFinder internal FASTAs to frozen original QfO inputs."""

import argparse
from collections import defaultdict
import json
from pathlib import Path
import re
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_input_sequences import sequence_identity
from benchmark_tools.inventory_swiss_sequences import PREPARED_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def read_maps(species_path, sequence_path):
    species, sequences, originals = {}, defaultdict(dict), set()
    for line in species_path.read_text().splitlines():
        key, separator, name = line.partition(": ")
        if not separator or not key.isdigit() or not name or Path(name).name != name or key in species or name in species.values():
            raise ValueError("Invalid or duplicate species map")
        species[key] = name
    with sequence_path.open() as stream:
        for line in stream:
            key, separator, description = line.rstrip("\n").partition(": ")
            match = re.fullmatch(r"(\d+)_(\d+)", key)
            if not separator or not match or match[1] not in species or not description.split():
                raise ValueError("Invalid sequence ID map")
            original = description.split()[0]
            bucket = sequences[match[1]]
            if key in bucket or original in originals:
                raise ValueError("Duplicate internal or original sequence ID")
            bucket[key] = original
            originals.add(original)
    if set(sequences) != set(species):
        raise ValueError("Empty or missing species sequence map")
    return species, sequences


def compare_species(original_path, internal_path, identifiers):
    original = {}
    for entry in SeqIO.parse(original_path, "fasta"):
        if entry.id in original:
            raise ValueError("Duplicate original FASTA ID")
        original[entry.id] = str(entry.seq)
    if set(original) != set(identifiers.values()) or len(original) != len(identifiers):
        raise ValueError("Original input identifiers differ from species map")
    seen, differences = set(), []
    for entry in SeqIO.parse(internal_path, "fasta"):
        if entry.id not in identifiers or entry.id in seen:
            raise ValueError("Unknown or duplicate internal FASTA ID")
        seen.add(entry.id)
        original_id = identifiers[entry.id]
        sequence = str(entry.seq)
        if sequence != original[original_id]:
            differences.append({"internal_id": entry.id, "original_id": original_id,
                                "original_sequence": sequence_identity(original[original_id]),
                                "internal_sequence": sequence_identity(sequence)})
    if seen != set(identifiers):
        raise ValueError("Incomplete internal sequence coverage")
    return {"sequences": len(original), "identical_sequences": len(original) - len(differences), "differences": differences}


def audit(prepared_path, working):
    prepared_identity = record(prepared_path)
    if prepared_identity["sha256"] != PREPARED_SHA:
        raise ValueError("Changed frozen original input manifest")
    prepared = json.loads(prepared_path.read_text())
    maps = [record(working / name) for name in ("SpeciesIDs.txt", "SequenceIDs.txt")]
    species, sequences = read_maps(*(Path(r["path"]) for r in maps))
    expected = {Path(r["path"]).name: r for r in prepared["input_fastas"]}
    if set(species.values()) != set(expected):
        raise ValueError("Proteome inventory differs")
    native_paths = {p.name for p in working.glob("Species*.fa")}
    if native_paths != {f"Species{index}.fa" for index in species}:
        raise ValueError("Unexpected internal FASTA inventory")
    rows, checked = [], [prepared_identity, *maps]
    for index, filename in species.items():
        native = record(working / f"Species{index}.fa")
        check(expected[filename])
        result = compare_species(Path(expected[filename]["path"]), Path(native["path"]), sequences[index])
        rows.append({"species_index": index, "original": expected[filename], "internal": native, **result})
        checked.extend((expected[filename], native))
    for identity in checked:
        check(identity)
    return {"status": "retained_orthofinder_input_sequences_compared", "source": record(__file__),
            "working_directory": str(working.resolve()), "checked_inputs": checked, "proteomes": rows,
            "total_sequences": sum(row["sequences"] for row in rows),
            "identical_sequences": sum(row["identical_sequences"] for row in rows),
            "sequence_difference_count": sum(len(row["differences"]) for row in rows),
            "limitations": ["Retained internal inputs checked, not an independently authenticated historical execution trace.",
                            "No inference, orthology scoring or proof that all other comparator inputs match.",
                            "This verifies parity with original frozen inputs, not corrected-release compatibility."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("prepared", "working", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.prepared, args.working)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")

"""Independently inventory staged corrected QfO FASTAs before inference freeze."""

import argparse
from bisect import bisect_right
import gzip
import json
from pathlib import Path
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.resolve_swiss_sequence_aliases import MAPPING_SHA


def inventory(directory, inputs, mapping):
    species, offsets = mapping["species"], mapping["Goff"]
    if (not species or len(set(species)) != len(species) or len(offsets) != len(species) + 1
            or offsets[0] != 0 or any(type(n) is not int for n in offsets)
            or any(a >= b for a, b in zip(offsets, offsets[1:]))):
        raise ValueError("Invalid reference species intervals")
    names = [Path(r["path"]).name for r in inputs]
    if len(names) != len(set(names)) or len(names) != len(species):
        raise ValueError("Incorrect staged proteome inventory")
    expected = set(names) | {"staging_manifest.json"}
    if {p.name for p in directory.iterdir()} != expected:
        raise ValueError("Unexpected staged directory entries")
    accessions, numeric_ids, owners, rows = set(), set(), set(), []
    for source in inputs:
        path = directory / Path(source["path"]).name
        if path.is_symlink() or not path.is_file() or str(path.resolve()) != source["path"]:
            raise ValueError("Staged FASTA is not a local regular file")
        check(source)
        counts, residues = {}, 0
        for entry in SeqIO.parse(path, "fasta"):
            parts = entry.id.split("|")
            sequence = str(entry.seq)
            if (len(parts) != 3 or parts[0] not in {"sp", "tr"} or not all(parts)
                    or not sequence or not sequence.isascii() or any(c.isspace() for c in sequence)):
                raise ValueError("Invalid staged sequence record")
            accession = parts[1]
            number = mapping["mapping"].get(accession)
            if accession in accessions or type(number) is not int or not 1 <= number <= offsets[-1] or number in numeric_ids:
                raise ValueError("Unmapped, duplicate or invalid staged identity")
            accessions.add(accession)
            numeric_ids.add(number)
            owner = species[bisect_right(offsets, number - 1) - 1]
            counts[owner] = counts.get(owner, 0) + 1
            residues += len(sequence)
        if len(counts) != 1:
            raise ValueError("Empty proteome or mixed reference species")
        owner = next(iter(counts))
        index = species.index(owner)
        if owner in owners or counts[owner] != offsets[index + 1] - offsets[index]:
            raise ValueError("Duplicate species or incomplete proteome coverage")
        owners.add(owner)
        rows.append({"file": source, "reference_species": owner,
                     "sequences": counts[owner], "residues": residues})
    if numeric_ids != set(range(1, offsets[-1] + 1)) or owners != set(species):
        raise ValueError("Incomplete reference universe")
    for source in inputs:
        path = directory / Path(source["path"]).name
        if path.is_symlink() or str(path.resolve()) != source["path"]:
            raise ValueError("Staged path changed during inventory")
        check(source)
    if {p.name for p in directory.iterdir()} != expected:
        raise ValueError("Staged directory changed during inventory")
    return {"proteomes": len(rows), "total_sequences": len(numeric_ids),
            "total_residues": sum(r["residues"] for r in rows), "files": rows,
            "unique_reference_species_per_proteome": True,
            "complete_unique_numeric_coverage": True}


def audit(manifest, manifest_sha, mapping_path):
    inputs = [record(manifest), record(mapping_path)]
    if [r["sha256"] for r in inputs] != [manifest_sha, MAPPING_SHA]:
        raise ValueError("Unreviewed staging manifest or changed reference mapping")
    stage = json.loads(manifest.read_text())
    if (manifest.name != "staging_manifest.json" or manifest.is_symlink()
            or stage["status"] != "corrected_inputs_staged_pending_inference_freeze"
            or stage["inputs_normalized"] is not False or stage["inference_authorized"] is not False
            or stage["total_sequences"] != 984137 or len(stage["input_fastas"]) != 78):
        raise ValueError("Unexpected corrected staging state")
    with gzip.open(mapping_path, "rt") as stream:
        mapping = json.load(stream)
    if mapping["Goff"][-1] != 984137 or len(mapping["species"]) != 78:
        raise ValueError("Unexpected reference universe")
    result = inventory(manifest.parent, stage["input_fastas"], mapping)
    for source in inputs:
        check(source)
    return {"status": "corrected_staged_inventory_verified_pending_execution_freeze",
            "source": record(__file__), "inputs": inputs, **result,
            "inference_authorized": False, "inputs_modified": False,
            "limitations": ["Reviewed staging-manifest hash binds prior compatibility review; this is not a substitute for it.",
                            "Direct FASTA inventory uses Biopython, also used by archive audits; not an independent parser implementation.",
                            "No inference, corrected accuracy, historical provenance or software-license approval established."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("manifest", "mapping", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.manifest, args.manifest_sha256, args.mapping)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(result["status"])

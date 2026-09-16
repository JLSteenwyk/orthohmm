"""Independently reconstruct the frozen YGOB reference without loading predictions."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import re
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance

# Zero-based columns transcribed independently from the acquired v7 README.
RETAINED = {0: "Vpolyspora", 1: "Tphaffii", 2: "Tblattae", 3: "Ndairenensis",
            4: "Ncastellii", 5: "Knaganishii", 6: "Kafricana", 7: "Cglabrata",
            13: "Zrouxii", 14: "Tdelbrueckii", 15: "Klactis", 16: "Egossypii",
            17: "Ecymbalariae", 18: "Lkluyveri", 19: "Lthermotolerans", 20: "Lwaltii",
            25: "Cglabrata", 26: "Kafricana", 27: "Knaganishii", 28: "Ncastellii",
            29: "Ndairenensis", 30: "Tblattae", 31: "Tphaffii", 32: "Vpolyspora"}


def reconstruct(pillars, fasta):
    with pillars.open() as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if any(len(row) != 33 for row in rows):
        raise ValueError("Expected 33 pillar columns")
    genes = [gene for row in rows for i, gene in enumerate(row) if i != 12 and gene != "---"]
    if any(not gene or any(c.isspace() for c in gene) for gene in genes):
        raise ValueError("Invalid pillar gene cell")
    duplicated = {g for g, n in Counter(genes).items() if n > 1}
    excluded_rows = {i for i, row in enumerate(rows, 1) if duplicated.intersection(row)}
    owners = {}
    for row in rows:
        for column, species in RETAINED.items():
            gene = row[column]
            if gene == "---":
                continue
            if gene in owners and owners[gene] != species:
                raise ValueError("Conflicting species assignment")
            owners[gene] = species
    sequences, seen = {}, set()
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id in seen:
            raise ValueError("Duplicate raw FASTA ID")
        seen.add(record.id)
        states = re.findall(r"\{(ON|OFF)\}", record.description)
        if len(states) != 1:
            raise ValueError("Missing or ambiguous ON/OFF state")
        if states[0] == "OFF" or record.id not in owners:
            continue
        sequence = str(record.seq).upper()
        if sequence.endswith("*"):
            sequence = sequence[:-1]
        if "*" in sequence:
            continue
        if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
            raise ValueError("Unsupported protein sequence")
        sequences[record.id] = (owners[record.id], sequence)
    reference = {}
    for index, row in enumerate(rows, 1):
        if index in excluded_rows:
            continue
        members = sorted(g for column, g in enumerate(row) if column in RETAINED and g in sequences)
        if members:
            reference[f"Pillar{index:05d}"] = members
    return sequences, reference, sorted(excluded_rows)


def verify(root):
    root = root.resolve()
    snapshot = root / "benchmark_tools/results/ygob_overlap_20260916.json"
    audit = json.loads(snapshot.read_text())
    for item in audit["candidate_inputs"]:
        verify_file(Path(item["path"]), item)
    candidate = root / "benchmarks/work/independent_ygob_v7"
    expected, reference, excluded = reconstruct(candidate / "Pillars.tab", candidate / "AA.fsa")
    prepared = root / "benchmarks/work/ygob_validation_v1"
    actual = {}
    files = sorted((prepared / "input").glob("*.fasta"))
    for path in files:
        for record in SeqIO.parse(path, "fasta"):
            if record.id in actual:
                raise ValueError("Duplicate prepared ID")
            actual[record.id] = (path.stem, str(record.seq))
    reference_path = prepared / "reference_groups.json"
    if expected != actual:
        raise ValueError("Prepared proteins/species/sequences differ from raw reconstruction")
    if reference != json.loads(reference_path.read_text()):
        raise ValueError("Prepared reference differs from raw reconstruction")
    reference_genes = {g for group in reference.values() for g in group}
    omitted = set(expected) - reference_genes
    if (len(actual), len(reference), len(reference_genes), len(omitted), excluded) != (83404, 10250, 83391, 13, [113, 9896]):
        raise ValueError("Frozen reference dimensions or exclusions changed")
    if {species for species, _ in actual.values()} != set(RETAINED.values()):
        raise ValueError("Frozen species set changed")
    return {"schema_version": 1, "status": "reference_reconstruction_verified", "accuracy_evaluated": False,
            "proteins": len(actual), "species": len(set(RETAINED.values())), "reference_groups": len(reference),
            "reference_genes": len(reference_genes), "excluded_rows": excluded,
            "input_genes_excluded_from_reference": sorted(omitted),
            "reference": file_provenance(reference_path), "source_inputs": audit["candidate_inputs"],
            "prepared_inputs": [file_provenance(p) for p in files], "verifier": file_provenance(Path(__file__)),
            "scope": "Independent raw-table and raw-FASTA reconstruction, not independent biological truth or prediction scoring."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = verify(args.root)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()

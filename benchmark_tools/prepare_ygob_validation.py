#!/usr/bin/env python3
"""Prepare label-independent YGOB group-recovery validation inputs."""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_ygob_overlap import PRE_WGD, POST_WGD, read_pillars
from benchmark_tools.orthobench_stage_diagnostics import file_provenance


EXCLUDED_SPECIES = ("Scerevisiae", "Skudriavzevii", "Smikatae", "Suvarum")
ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def prepare(candidate, output, audit_path):
    snapshot = json.loads(audit_path.read_text())
    for source in snapshot["candidate_inputs"]:
        actual = file_provenance(candidate / Path(source["path"]).name)
        if actual["sha256"] != source["sha256"]:
            raise ValueError(f"Candidate snapshot changed: {actual['path']}")
    _, assignments = read_pillars(candidate / "Pillars.tab")
    ambiguous_rows = {line for r in assignments.values() if len(r["pillar_lines"]) > 1 for line in r["pillar_lines"]}
    excluded_reference = {g for g, r in assignments.items() if ambiguous_rows.intersection(r["pillar_lines"])}
    records = defaultdict(list)
    references = defaultdict(list)
    exclusions = {"OFF": [], "exposed_genus": [], "internal_stop": []}
    observed = set()
    for record in SeqIO.parse(candidate / "AA.fsa", "fasta"):
        if record.id in observed:
            raise ValueError(f"Duplicate protein: {record.id}")
        observed.add(record.id)
        states = re.findall(r"\{(ON|OFF)\}", record.description)
        if len(states) != 1:
            raise ValueError(f"Invalid state annotation: {record.id}")
        if states[0] == "OFF":
            exclusions["OFF"].append(record.id)
            continue
        if record.id not in assignments:
            raise ValueError(f"Unmapped ON protein: {record.id}")
        assignment = assignments[record.id]
        species = assignment["species"]
        if species in EXCLUDED_SPECIES:
            exclusions["exposed_genus"].append(record.id)
            continue
        sequence = str(record.seq).upper().removesuffix("*")
        if "*" in sequence:
            exclusions["internal_stop"].append(record.id)
            continue
        if not sequence or set(sequence) - ALPHABET:
            raise ValueError(f"Invalid protein sequence: {record.id}")
        records[species].append(SeqRecord(Seq(sequence), id=record.id, description=""))
        if record.id not in excluded_reference:
            references[f"Pillar{assignment['pillar_lines'][0]:05d}"].append(record.id)
    if not records:
        raise ValueError("No retained proteins")
    if output.exists():
        raise ValueError(f"Refusing to overwrite prepared inputs: {output}")
    output.mkdir(parents=True)
    inputs = output / "input"
    inputs.mkdir()
    for species, values in sorted(records.items()):
        SeqIO.write(sorted(values, key=lambda r: r.id), inputs / (species + ".fasta"), "fasta")
    reference_path = output / "reference_groups.json"
    reference_path.write_text(json.dumps({k: sorted(v) for k, v in sorted(references.items())}, indent=2) + "\n")
    input_genes = {r.id for values in records.values() for r in values}
    result = {
        "schema_version": 1, "generated_at": datetime.now(timezone.utc).isoformat(),
        "description": "YGOB v7 group-recovery candidate; no accuracy results",
        "source": file_provenance(Path(__file__)), "candidate_audit": file_provenance(audit_path),
        "excluded_species": list(EXCLUDED_SPECIES),
        "species": sorted(records), "proteins_by_species": {s: len(v) for s, v in sorted(records.items())},
        "proteins": len(input_genes), "reference_groups": len(references),
        "reference_genes": sum(map(len, references.values())),
        "excluded_ambiguous_pillar_rows": sorted(ambiguous_rows),
        "input_genes_excluded_from_reference": sorted(input_genes & excluded_reference),
        "input_exclusions": {k: sorted(v) for k, v in exclusions.items()},
        "inputs": [file_provenance(p) for p in sorted(inputs.glob("*.fasta"))],
        "reference": file_provenance(reference_path),
        "interpretation": "Curated homology group recovery; not A/B-resolved post-WGD pairwise orthology",
        "independent_family_validation_established": False,
    }
    (output / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.candidate, args.output, args.audit)
    expected = set(PRE_WGD + POST_WGD) - set(EXCLUDED_SPECIES)
    if set(result["species"]) != expected:
        raise ValueError("Prepared species do not match the prespecified 16-species set")
    result["command"] = [sys.executable, *sys.argv]
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"Prepared {len(result['species'])} species, {result['proteins']} proteins, {result['reference_groups']} reference groups")


if __name__ == "__main__":
    main()

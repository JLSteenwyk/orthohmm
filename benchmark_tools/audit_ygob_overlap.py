#!/usr/bin/env python3
"""Inventory YGOB v7 proteins and exact overlap before any accuracy evaluation."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.orthobench_stage_diagnostics import file_provenance


# Documented v7 pillar columns, not inferred A/B subgenome assignments.
POST_WGD = ("Vpolyspora", "Tphaffii", "Tblattae", "Ndairenensis", "Ncastellii",
            "Knaganishii", "Kafricana", "Cglabrata", "Suvarum", "Skudriavzevii",
            "Smikatae", "Scerevisiae")
PRE_WGD = ("Zrouxii", "Tdelbrueckii", "Klactis", "Egossypii", "Ecymbalariae",
           "Lkluyveri", "Lthermotolerans", "Lwaltii")
COLUMNS = POST_WGD + ("Ancestor",) + PRE_WGD + POST_WGD[::-1]


def fingerprint(sequence):
    # Ignore case and one terminal translation stop only; preserve internal stops.
    return hashlib.sha256(str(sequence).upper().removesuffix("*").encode("ascii")).hexdigest()


def read_pillars(path):
    genes = {}
    count = 0
    with path.open() as handle:
        for count, fields in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(fields) != len(COLUMNS):
                raise ValueError(f"Expected 33 columns at {path}:{count}, got {len(fields)}")
            for species, gene in zip(COLUMNS, fields):
                if gene == "---" or species == "Ancestor":
                    continue
                if not gene or gene.strip() != gene or any(c.isspace() for c in gene):
                    raise ValueError(f"Invalid gene cell at {path}:{count}")
                if gene in genes:
                    if genes[gene]["species"] != species:
                        raise ValueError(f"Conflicting species for pillar gene: {gene}")
                    genes[gene]["pillar_lines"].append(count)
                else:
                    genes[gene] = {"species": species, "pillar_lines": [count]}
    return count, genes


def audit(candidate, development):
    pillar_count, genes = read_pillars(candidate / "Pillars.tab")
    sequences = {}
    by_hash = defaultdict(list)
    off = 0
    seen = set()
    missing = []
    for record in SeqIO.parse(candidate / "AA.fsa", "fasta"):
        if record.id in seen:
            raise ValueError(f"Duplicate protein ID: {record.id}")
        seen.add(record.id)
        states = re.findall(r"\{(ON|OFF)\}", record.description)
        if len(states) != 1:
            raise ValueError(f"Missing or ambiguous ON/OFF annotation: {record.id}")
        if states[0] == "OFF":
            off += 1
            continue
        if record.id not in genes:
            missing.append(record.id)
        sequence = str(record.seq).upper().removesuffix("*")
        if not sequence:
            raise ValueError(f"Empty protein: {record.id}")
        digest = fingerprint(record.seq)
        assignment = genes.get(record.id, {"species": "unassigned", "pillar_lines": []})
        sequences[record.id] = {**assignment, "sha256": digest, "length": len(sequence),
                                "internal_stop": "*" in sequence}
        by_hash[digest].append(record.id)
    if not sequences:
        raise ValueError("No ON proteins found")
    overlaps = {}
    for label, directory in development.items():
        files = sorted(p for p in directory.iterdir() if p.suffix in {".fa", ".faa", ".fasta", ".fsa"})
        if not files:
            raise ValueError(f"No development FASTAs: {directory}")
        matched = set()
        total = matching_records = 0
        for path in files:
            for record in SeqIO.parse(path, "fasta"):
                total += 1
                matching = by_hash.get(fingerprint(record.seq), ())
                if matching:
                    matching_records += 1
                    matched.update(matching)
        overlaps[label] = {
            "development_proteins": total, "matching_development_records": matching_records,
            "matching_candidate_proteins": len(matched),
            "matching_candidate_pillars": len({p for g in matched for p in sequences[g]["pillar_lines"]}),
            "candidate_matches_by_species": dict(Counter(sequences[g]["species"] for g in matched)),
            "matching_candidate_gene_ids": sorted(matched),
            "inputs": [file_provenance(p) for p in files],
        }
    return {
        "schema_version": 1, "candidate": "YGOB v7-Aug2012", "scoring_performed": False,
        "independence_established": False, "pillar_rows": pillar_count,
        "pillar_gene_slots": sum(len(r["pillar_lines"]) for r in genes.values()),
        "unique_pillar_gene_ids": len(genes), "on_proteins": len(sequences), "off_proteins": off,
        "on_proteins_without_pillar": sorted(missing),
        "repeated_pillar_gene_ids": {g: r for g, r in genes.items() if len(r["pillar_lines"]) > 1},
        "on_proteins_by_species": dict(Counter(r["species"] for r in sequences.values())),
        "internal_stop_proteins": sorted(g for g, r in sequences.items() if r["internal_stop"]),
        "pillar_genes_without_protein_sequence": len(set(genes) - set(sequences)),
        "overlaps": overlaps,
        "limitations": [
            "Exact sequence matches establish overlap; absence of exact matches does not establish family independence.",
            "Species identity and synonyms still require taxonomic auditing; tiny exact matches can occur across species.",
            "YGOB pillars are curated homology groups, not A/B-resolved post-WGD pairwise orthology labels.",
            "Non-protein features can occur in pillars; missing protein sequences are not automatically annotation failures.",
            "Repeated pillar membership and unmapped ON proteins must be resolved before using a partition-based scorer.",
            "Data were retrieved from the official HTTP archive; checksums identify the snapshot but do not authenticate transport.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--development", action="append", required=True, help="NAME=DIRECTORY")
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    development = {}
    for entry in args.development:
        label, path = entry.split("=", 1)
        if not label or label in development:
            raise ValueError("Development names must be nonempty and unique")
        development[label] = Path(path)
    result = audit(args.candidate, development)
    result.update(generated_at=datetime.now(timezone.utc).isoformat(), command=[sys.executable, *sys.argv],
                  source=file_provenance(Path(__file__)),
                  candidate_inputs=[file_provenance(args.candidate / n) for n in ("AA.fsa", "Pillars.tab", "README")])
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"{result['on_proteins']} ON proteins; {result['pillar_rows']} pillar rows")
    for label, values in result["overlaps"].items():
        print(f"{label}: {values['matching_candidate_proteins']} candidate proteins match exactly")


if __name__ == "__main__":
    main()

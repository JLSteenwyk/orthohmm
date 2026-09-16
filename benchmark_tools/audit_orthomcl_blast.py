#!/usr/bin/env python3
"""Audit legacy BLAST diagnostics against inputs and final OrthoMCL groups."""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.qfo_filter_pairs import load_mapping


DIAGNOSTIC = re.compile(
    r"^\[blastall\] (WARNING|ERROR):\s+\[[^\]]+\]\s+(\S+): (.+)$"
)


def parse_diagnostics(path):
    records = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        match = DIAGNOSTIC.fullmatch(line)
        if not match:
            raise ValueError(f"Unrecognized BLAST diagnostic at {path}:{number}")
        level, gene, message = match.groups()
        if message == "SetUpBlastSearch failed.":
            category = "setup_failure"
        elif "Unable to calculate Karlin-Altschul params" in message:
            category = "statistics_failure"
        elif "Query must be at least twice wordsize" in message:
            category = "short_query_failure"
        elif "Selenocysteine (U)" in message and "replaced by X" in message:
            category = "selenocysteine_replacement"
        else:
            raise ValueError(f"Unclassified BLAST diagnostic at {path}:{number}")
        record = records.setdefault(gene, {"gene": gene, "messages": []})
        record["messages"].append({
            "line": number, "level": level, "category": category, "message": message,
        })
    for record in records.values():
        record["query_failed"] = any(
            m["category"].endswith("failure") for m in record["messages"]
        )
    return records


def audit(log, input_dir, groups, mapping):
    records = parse_diagnostics(log)
    seen = set()
    species_counts = Counter()
    fasta_paths = sorted(input_dir.glob("*.fasta"))
    if not fasta_paths:
        raise ValueError(f"No FASTA inputs in {input_dir}")
    for path in fasta_paths:
        for seq in SeqIO.parse(path, "fasta"):
            if seq.id in seen:
                raise ValueError(f"Duplicate FASTA ID: {seq.id}")
            seen.add(seq.id)
            species_counts[path.stem] += 1
            if seq.id in records:
                record = records[seq.id]
                record.update(species=path.stem, length=len(seq))
                record["residue_counts"] = dict(Counter(str(seq.seq).upper()))
                accession = seq.id.split("|")[1] if "|" in seq.id else seq.id
                record["accession"] = accession
                record["qfo_mapping_valid"] = accession in mapping
                record["final_group_line"] = None
                record["final_group_size"] = 0
    missing = set(records) - seen
    if missing:
        raise ValueError(f"Diagnostic IDs absent from inputs: {sorted(missing)}")
    grouped = set()
    group_count = 0
    with groups.open() as handle:
        for number, line in enumerate(handle, 1):
            genes = line.split()
            if not genes:
                continue
            if len(set(genes)) != len(genes) or grouped.intersection(genes):
                raise ValueError(f"Duplicate group membership at {groups}:{number}")
            if not set(genes) <= seen:
                raise ValueError(f"Unknown group member at {groups}:{number}")
            grouped.update(genes)
            group_count += 1
            for gene in set(genes).intersection(records):
                records[gene]["final_group_line"] = number
                records[gene]["final_group_size"] = len(genes)
    failed = [r for r in records.values() if r["query_failed"]]
    replaced = [r for r in records.values() if any(
        m["category"] == "selenocysteine_replacement" for m in r["messages"]
    )]
    return {
        "schema_version": 1,
        "input_proteins": len(seen),
        "input_species": len(fasta_paths),
        "input_proteins_by_species": dict(species_counts),
        "final_groups": group_count,
        "final_grouped_proteins": len(grouped),
        "diagnostic_lines_by_category": dict(Counter(
            m["category"] for r in records.values() for m in r["messages"]
        )),
        "failed_queries": len(failed),
        "failed_query_fraction": len(failed) / len(seen),
        "failed_queries_in_final_groups": sum(r["final_group_line"] is not None for r in failed),
        "failed_queries_qfo_mapping_valid": sum(r["qfo_mapping_valid"] for r in failed),
        "failed_queries_by_species": dict(Counter(r["species"] for r in failed)),
        "selenocysteine_replaced_queries": len(replaced),
        "selenocysteine_replaced_queries_in_final_groups": sum(
            r["final_group_line"] is not None for r in replaced
        ),
        "records": [records[gene] for gene in sorted(records)],
        "limitations": [
            "Query setup failures do not prove absence of incoming subject hits.",
            "Mapping membership does not establish membership in a QfO reference family.",
            "Final group membership does not establish correctness of the assignment.",
            "This audit does not yet measure the counterfactual score after repairing failed queries.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--groups", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    args = parser.parse_args()
    result = audit(args.log, args.input_dir, args.groups, load_mapping(args.mapping))
    result.update(
        generated_at=datetime.now(timezone.utc).isoformat(),
        command=[sys.executable, *sys.argv],
        source=file_provenance(Path(__file__)),
        inputs={
            "log": file_provenance(args.log),
            "groups": file_provenance(args.groups),
            "mapping": file_provenance(args.mapping),
            "proteomes": [file_provenance(p) for p in sorted(args.input_dir.glob("*.fasta"))],
        },
    )
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in {
        "records", "inputs", "input_proteins_by_species", "source", "command"
    }}, indent=2))


if __name__ == "__main__":
    main()

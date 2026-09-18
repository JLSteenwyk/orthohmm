"""Locate missing SwissTrees accessions in retained canonical/additional FASTAs."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.resolve_swiss_sequence_aliases import INVENTORY_SHA, MAPPING_SHA
from benchmark_tools.qfo_filter_pairs import load_mapping

ALIAS_SHA = "b69f74e35aa3fa3af99a33b3b27f830d79df2abd3935b84e911f981186ac33ec"


def scan(path, selected=None):
    found = {}
    count = 0
    for entry in SeqIO.parse(path, "fasta"):
        parts = entry.id.split("|")
        if len(parts) != 3 or parts[0] not in {"sp", "tr"} or not all(parts):
            raise ValueError("Unexpected FASTA header")
        accession = parts[1]
        count += 1
        if selected is not None and accession not in selected:
            continue
        if accession in found:
            raise ValueError("Duplicate selected accession")
        match = re.match(r"\S+ Isoform of ([A-Za-z0-9]+),", entry.description)
        sequence = str(entry.seq)
        if not sequence:
            raise ValueError("Empty selected sequence")
        found[accession] = {"description": entry.description, "length": len(sequence),
                            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                            "header_isoform_of": match.group(1) if match else None}
    return count, found


def audit(alias_path, canonical, additional, staged, mapping_path):
    sources = [record(p) for p in (alias_path, canonical, additional, staged, mapping_path)]
    if sources[0]["sha256"] != ALIAS_SHA or sources[-1]["sha256"] != MAPPING_SHA:
        raise ValueError("Changed frozen audit/mapping")
    aliases = json.loads(alias_path.read_text())
    staged_records = [r for r in aliases["fasta_inputs"] if r["path"] == str(staged.resolve())]
    if staged_records != [sources[3]]:
        raise ValueError("Staged FASTA is not the frozen input")
    if (sources[1]["bytes"], sources[1]["sha256"]) != (sources[3]["bytes"], sources[3]["sha256"]):
        raise ValueError("Canonical and staged FASTAs differ")
    missing = set(aliases["summary"]["missing_genes"])
    canonical_count, primary = scan(canonical)
    additional_count, extra = scan(additional, missing)
    if missing & primary.keys():
        raise ValueError("Missing identity appears in canonical input")
    mapping = load_mapping(mapping_path)
    rows = {}
    for gene in sorted(missing):
        value = extra.get(gene)
        target = value["header_isoform_of"] if value else None
        rows[gene] = {"found_in_additional": value is not None, "additional_record": value,
                      "numeric_protein_id": aliases["identities"][gene]["numeric_protein_id"],
                      "header_target_in_canonical": target in primary if target else None,
                      "header_target_numeric_id": mapping.get(target) if target else None,
                      "header_target_record": primary.get(target) if target else None}
    for identity in sources:
        check(identity)
    return {"status": "retained_additional_sequence_inventory", "source": record(__file__),
            "inputs": sources, "canonical_staged_byte_identical": True,
            "canonical_sequences": canonical_count, "additional_sequences": additional_count,
            "missing_reference_genes": len(missing), "found_in_additional": len(extra), "genes": rows,
            "limitations": ["Retained extracted files checked; archive-member authentication is not performed here.",
                            "Isoform-of text is source annotation, not proven orthology or authorization to remap reference truth.",
                            "No sequence-similarity search, scoring, input changes or historical comparator input reconstruction performed."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("aliases", "canonical", "additional", "staged", "mapping", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.aliases, args.canonical, args.additional, args.staged, args.mapping)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")

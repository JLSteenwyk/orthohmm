"""Resolve SwissTrees sequence identities through the frozen QfO numeric mapping."""

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.inventory_swiss_sequences import collect, summarize, PREPARED_SHA
from benchmark_tools.qfo_filter_pairs import load_mapping
from benchmark_tools.snapshot_orthohmm_input_order import record

INVENTORY_SHA = "f4a147ff1ecdf0c8df046ab1a5c63069fef106fa627d52eaec008aae250527a5"
MAPPING_SHA = "1c10f6ce5e53ebc3148dde02d16268b225c3c8817c952d1274daa41acbf9eb4d"


def resolve(references, input_accessions, mapping):
    wanted = {}
    for accession in references:
        number = mapping.get(accession)
        if type(number) is not int or number <= 0:
            raise ValueError("Reference accession lacks valid numeric identity")
        if number in wanted:
            raise ValueError("Shared reference numeric identity")
        wanted[number] = accession
    candidates = defaultdict(list)
    seen = set()
    for accession in input_accessions:
        if accession in seen:
            raise ValueError("Duplicate input accession")
        seen.add(accession)
        number = mapping.get(accession)
        if number is not None and (type(number) is not int or number <= 0):
            raise ValueError("Invalid input numeric identity")
        if number in wanted:
            candidates[number].append(accession)
    rows = {}
    for number, accession in sorted(wanted.items()):
        matches = sorted(candidates[number])
        status = "missing" if not matches else "ambiguous" if len(matches) > 1 else "exact" if matches[0] == accession else "mapped_alias"
        rows[accession] = {"numeric_protein_id": number, "input_accessions": matches, "status": status}
    return rows


def audit(inventory_path, prepared_path, mapping_path):
    sources = [record(p) for p in (inventory_path, prepared_path, mapping_path)]
    if [r["sha256"] for r in sources] != [INVENTORY_SHA, PREPARED_SHA, MAPPING_SHA]:
        raise ValueError("Changed frozen identity source")
    inventory = json.loads(inventory_path.read_text())
    prepared = json.loads(prepared_path.read_text())
    if inventory["fasta_inputs"] != prepared["input_fastas"]:
        raise ValueError("Different input inventories")
    count_path = Path(inventory["inputs"][0]["path"])
    if record(count_path) != inventory["inputs"][0]:
        raise ValueError("Changed reference membership")
    counts = json.loads(count_path.read_text())
    families = {r["family"]: r["represented_genes"] for r in counts["methods"][0]["families"]}
    accessions = []
    for identity in prepared["input_fastas"]:
        path = Path(identity["path"])
        if record(path) != identity:
            raise ValueError("Changed FASTA input")
        for entry in SeqIO.parse(path, "fasta"):
            parts = entry.id.split("|")
            if len(parts) != 3 or parts[0] not in {"sp", "tr"} or not all(parts):
                raise ValueError("Unexpected FASTA identifier")
            accessions.append(parts[1])
        if record(path) != identity:
            raise ValueError("Input changed during identity indexing")
    references = [g for members in families.values() for g in members]
    identities = resolve(references, accessions, load_mapping(mapping_path))
    selected = {g: r["input_accessions"][0] for g, r in identities.items() if r["status"] in {"exact", "mapped_alias"}}
    mapped = collect({f: [selected[g] for g in members if g in selected] for f, members in families.items()}, prepared["input_fastas"])
    found = {g: {**mapped["genes"][a], "resolved_input_accession": a, "identity_resolution": identities[g]} for g, a in selected.items()}
    for gene, old in inventory["genes"].items():
        if identities[gene]["status"] != "exact" or any(found[gene][key] != value for key, value in old.items()):
            raise ValueError("Original exact-match descriptor changed")
    for identity in [*sources, inventory["inputs"][0]]:
        if record(identity["path"]) != identity:
            raise ValueError("Source changed during audit")
    return {"status": "swiss_sequence_numeric_identity_audit", "source": record(__file__),
            "helpers": [record(Path(__file__).with_name(name)) for name in
                        ("inventory_swiss_sequences.py", "qfo_filter_pairs.py", "snapshot_orthohmm_input_order.py")],
            "inputs": sources, "fasta_inputs": prepared["input_fastas"],
            "resolution_counts": {status: sum(r["status"] == status for r in identities.values())
                                  for status in ("exact", "mapped_alias", "missing", "ambiguous")},
            "identities": identities, "genes": found, "summary": summarize(references, found),
            "families": {f: summarize(members, found) for f, members in families.items()},
            "prediction_statistics_evaluated": False,
            "limitations": ["Numeric identity is defined by the frozen QfO mapping, not sequence similarity or current external annotation.",
                            "Mapped aliases identify benchmark inputs; this does not prove historical alias sequences were byte-identical.",
                            "Missing and ambiguous identities are unresolved, never imputed.",
                            "Original entropy and fragment-label limitations remain; no absence-of-fragments or causal claim."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("inventory", "prepared", "mapping", "output"):
        parser.add_argument("--" + flag, required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.inventory, args.prepared, args.mapping)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")

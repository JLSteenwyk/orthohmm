"""Inventory SwissTrees input sequence descriptors without evaluating predictions."""

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import statistics
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.inventory_swiss_annotations import COUNTS_SHA
from benchmark_tools.snapshot_orthohmm_input_order import record

PREPARED_SHA = "706b07c91e9a130dae229837641a7ad62d7d36a09679e5e0daa0959e182b7d64"
AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def describe(sequence, description):
    sequence = sequence.upper()
    if not sequence or any(c.isspace() for c in sequence):
        raise ValueError("Empty or whitespace-containing sequence")
    counts = Counter(sequence)
    canonical = sum(counts[c] for c in AMINO_ACIDS)
    frequencies = [counts[c] / canonical for c in sorted(AMINO_ACIDS) if counts[c]] if canonical else []
    return {
        "length": len(sequence),
        "canonical_residues": canonical,
        "noncanonical_counts": {c: n for c, n in sorted(counts.items()) if c not in AMINO_ACIDS},
        "noncanonical_fraction": (len(sequence) - canonical) / len(sequence),
        "canonical_entropy_bits": -sum(p * math.log2(p) for p in frequencies) if canonical else None,
        "maximum_canonical_frequency": max(frequencies) if canonical else None,
        "explicit_fragment_description": any(marker in description.lower() for marker in ("(fragment)", "(fragments)")),
        "description": description,
    }


def summarize(members, found):
    selected = [found[g] for g in members if g in found]
    entropy = [r["canonical_entropy_bits"] for r in selected if r["canonical_entropy_bits"] is not None]
    return {
        "reference_genes": len(members), "matched_genes": len(selected),
        "missing_genes": sorted(set(members) - set(found)),
        "explicit_fragment_descriptions": sum(r["explicit_fragment_description"] for r in selected),
        "genes_with_noncanonical_residues": sum(r["noncanonical_fraction"] > 0 for r in selected),
        "genes_without_canonical_residues": sum(r["canonical_residues"] == 0 for r in selected),
        "median_canonical_entropy_bits": statistics.median(entropy) if entropy else None,
        "median_length": statistics.median(r["length"] for r in selected) if selected else None,
    }


def collect(families, inputs):
    genes = [g for members in families.values() for g in members]
    if len(genes) != len(set(genes)):
        raise ValueError("Shared or duplicate reference gene")
    wanted, found = set(genes), {}
    if len({r["path"] for r in inputs}) != len(inputs):
        raise ValueError("Duplicate input path")
    for identity in inputs:
        path = Path(identity["path"])
        if record(path) != identity:
            raise ValueError("Changed FASTA input")
        for entry in SeqIO.parse(path, "fasta"):
            fields = entry.id.split("|")
            if len(fields) != 3 or fields[0] not in {"sp", "tr"} or not all(fields):
                raise ValueError("Unexpected accession header")
            accession = fields[1]
            if accession not in wanted:
                continue
            if accession in found:
                raise ValueError("Ambiguous reference accession")
            found[accession] = {"source_file": path.name, "input_id": entry.id,
                                **describe(str(entry.seq), entry.description)}
        if record(path) != identity:
            raise ValueError("FASTA changed during inventory")
    return {"summary": summarize(genes, found), "genes": found,
            "families": {f: summarize(members, found) for f, members in families.items()}}


def inventory(count_path, prepared_path):
    sources = [record(count_path), record(prepared_path)]
    if [r["sha256"] for r in sources] != [COUNTS_SHA, PREPARED_SHA]:
        raise ValueError("Changed frozen membership or input manifest")
    counts, prepared = json.loads(count_path.read_text()), json.loads(prepared_path.read_text())
    families = {r["family"]: r["represented_genes"] for r in counts["methods"][0]["families"]}
    if len(families) != 18 or sum(map(len, families.values())) != 563:
        raise ValueError("Unexpected SwissTrees reference universe")
    result = collect(families, prepared["input_fastas"])
    if [record(count_path), record(prepared_path)] != sources:
        raise ValueError("Manifest changed during inventory")
    return {"status": "prediction_independent_swiss_sequence_inventory", **result,
            "source": record(__file__), "inputs": sources, "fasta_inputs": prepared["input_fastas"],
            "prediction_statistics_evaluated": False,
            "limitations": [
                "Composition is global canonical-residue Shannon entropy, not local low-complexity or a causal mechanism.",
                "Noncanonical symbols are reported separately and excluded from the entropy denominator.",
                "Fragment labels mean literal parenthesized Fragment/Fragments in retained FASTA descriptions only.",
                "An absent fragment label does not establish complete sequence, and labels are not independently validated truncations.",
                "No score-based strata, thresholds, contrasts or uncertainty intervals evaluated."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("counts", "prepared", "output"):
        parser.add_argument("--" + flag, required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = inventory(args.counts, args.prepared)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")

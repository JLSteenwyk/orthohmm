"""Compare retained OrthoMCL combined inputs to frozen original QfO FASTAs."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_input_sequences import sequence_identity
from benchmark_tools.inventory_swiss_sequences import PREPARED_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def compare_inputs(sources, combined, genome_map):
    expected, taxa = {}, set()
    for source in sources:
        taxon = source.stem
        if taxon in taxa:
            raise ValueError("Duplicate original taxon")
        taxa.add(taxon)
        count = 0
        for entry in SeqIO.parse(source, "fasta"):
            if entry.id in expected:
                raise ValueError("Duplicate original FASTA ID")
            if not entry.seq:
                raise ValueError("Empty original sequence")
            expected[entry.id] = (taxon, sequence_identity(str(entry.seq)))
            count += 1
        if not count:
            raise ValueError("Empty original proteome")
    seen, differences, counts = set(), [], Counter()
    for entry in SeqIO.parse(combined, "fasta"):
        if entry.id not in expected or entry.id in seen:
            raise ValueError("Unknown or duplicate combined FASTA ID")
        seen.add(entry.id)
        taxon, original = expected[entry.id]
        counts[taxon] += 1
        native = sequence_identity(str(entry.seq))
        if native != original:
            differences.append({"id": entry.id, "taxon": taxon,
                                "original": original, "combined": native})
    if seen != set(expected):
        raise ValueError("Incomplete combined sequence coverage")
    mapped, mapped_taxa = set(), set()
    with genome_map.open() as stream:
        for line in stream:
            taxon, separator, genes = line.rstrip("\n").partition(":")
            identifiers = genes.split()
            if not separator or taxon not in taxa or taxon in mapped_taxa or not identifiers:
                raise ValueError("Invalid or duplicate genome map taxon")
            mapped_taxa.add(taxon)
            for identifier in identifiers:
                if identifier not in expected or identifier in mapped:
                    raise ValueError("Unknown or duplicate genome map ID")
                if expected[identifier][0] != taxon:
                    raise ValueError("Incorrect genome map species assignment")
                mapped.add(identifier)
    if mapped != set(expected) or mapped_taxa != taxa:
        raise ValueError("Incomplete genome map coverage")
    return {"proteomes": len(taxa), "total_sequences": len(expected),
            "identical_sequences": len(expected) - len(differences),
            "sequence_difference_count": len(differences), "differences": differences,
            "sequences_per_taxon": dict(sorted(counts.items())),
            "genome_map_complete_and_correct": True}


def audit(prepared, combined, genome_map):
    identity = record(prepared)
    if identity["sha256"] != PREPARED_SHA:
        raise ValueError("Changed frozen original input manifest")
    inputs = json.loads(prepared.read_text())["input_fastas"]
    checked = [identity, *inputs, record(combined), record(genome_map)]
    for item in checked:
        check(item)
    result = compare_inputs([Path(item["path"]) for item in inputs], combined, genome_map)
    for item in checked:
        check(item)
    return {"status": "retained_orthomcl_input_sequences_compared", "source": record(__file__),
            "checked_inputs": checked, **result,
            "limitations": ["Retained FASTA and genome map, not an authenticated historical execution trace or BLAST database audit.",
                            "Exact case-sensitive sequence length and SHA256 comparison; no residue normalization.",
                            "Parity with original frozen inputs does not establish corrected-release compatibility.",
                            "Does not resolve the documented sequence-specific BLAST failures or establish other tools' input parity."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("prepared", "combined", "genome-map", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.prepared, args.combined, args.genome_map)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")

"""Prepare complete four-species WGD inputs and outcome-blind reference mapping."""

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_ygob_overlap import read_pillars
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_ygob_validation import ALPHABET
from benchmark_tools.snapshot_orthohmm_input_order import record

SPECIES = ("Scerevisiae", "Skudriavzevii", "Smikatae", "Suvarum")


def collect_sequences(path, assignments):
    records = {species: [] for species in SPECIES}
    status, excluded = {}, defaultdict(list)
    for entry in SeqIO.parse(path, "fasta"):
        if entry.id in status:
            raise ValueError("Duplicate source protein ID")
        states = re.findall(r"\{(ON|OFF)\}", entry.description)
        if len(states) != 1:
            raise ValueError("Missing or ambiguous source ON/OFF state")
        if states[0] == "OFF":
            status[entry.id] = "OFF"
            excluded["OFF"].append(entry.id)
            continue
        if entry.id not in assignments:
            raise ValueError("ON protein lacks source species assignment")
        species = assignments[entry.id]["species"]
        if species not in SPECIES:
            status[entry.id] = "outside_four_species"
            continue
        sequence = str(entry.seq).upper().removesuffix("*")
        if "*" in sequence:
            status[entry.id] = "internal_stop"
            excluded["internal_stop"].append(entry.id)
            continue
        if not sequence or set(sequence) - ALPHABET:
            raise ValueError("Invalid retained protein sequence")
        status[entry.id] = "included"
        records[species].append(SeqRecord(Seq(sequence), id=entry.id, description=""))
    if any(not rows for rows in records.values()):
        raise ValueError("Missing required complete-proteome species")
    return records, status, dict(excluded)


def map_cohort(cohort, assignments, owners, status):
    ambiguous = {line for value in assignments.values() if len(value["pillar_lines"]) != 1
                 for line in value["pillar_lines"]}
    reference = defaultdict(list)
    for gene in owners:
        value = assignments[gene]
        if len(value["pillar_lines"]) == 1 and value["pillar_lines"][0] not in ambiguous:
            reference[f"Pillar{value['pillar_lines'][0]:05d}"].append(gene)
    mapped = []
    for pair in cohort:
        input_reasons, reference_reasons, lines = [], [], []
        for gene in pair["orf_pair"]:
            if gene not in owners:
                input_reasons.append({"gene": gene, "reason": status.get(gene, "no_protein_sequence")})
            elif owners[gene] != "Scerevisiae":
                input_reasons.append({"gene": gene, "reason": "unexpected_species"})
            value = assignments.get(gene)
            if value is None:
                reference_reasons.append({"gene": gene, "reason": "no_pillar_assignment"})
            elif len(value["pillar_lines"]) != 1 or value["pillar_lines"][0] in ambiguous:
                reference_reasons.append({"gene": gene, "reason": "ambiguous_pillar"})
            else:
                lines.append(value["pillar_lines"][0])
        if len(lines) == 2 and lines[0] != lines[1]:
            reference_reasons.append({"reason": "anchors_in_different_pillars", "lines": lines})
        pillar = f"Pillar{lines[0]:05d}" if len(lines) == 2 and lines[0] == lines[1] and not reference_reasons else None
        members = sorted(reference[pillar]) if pillar is not None else []
        mapped.append({**pair, "split_eligible": not input_reasons,
                       "reference_eligible": not input_reasons and not reference_reasons,
                       "input_reasons": input_reasons, "reference_reasons": reference_reasons,
                       "reference_pillar": pillar, "available_pillar_members": members,
                       "available_members_by_species": dict(Counter(owners[g] for g in members))})
    return mapped, {key: sorted(value) for key, value in sorted(reference.items())}, sorted(ambiguous)


def prepare(candidate, cohort_path, audit_path, protocol, output):
    cohort = read_pinned(cohort_path, "46f0bdf5bae23c8d30660aac59ca1a71a26deec5753a8d24ba09e04f86919ac5")
    audit = read_pinned(audit_path, "7508f53326bd884d97bca92e0bcb49a56e252e7a0c9333361ac46d820fed9204")
    if record(protocol)["sha256"] != "8126f75eaf30c34d4988c73ea233127a1f7f069fcb8487abaafaa1cf5200569d":
        raise ValueError("Changed application protocol")
    sources = []
    for expected in audit["candidate_inputs"]:
        observed = record(candidate / Path(expected["path"]).name)
        if observed["sha256"] != expected["sha256"]:
            raise ValueError("Changed original YGOB snapshot")
        sources.append(observed)
    _, assignments = read_pillars(candidate / "Pillars.tab")
    sequences, status, excluded = collect_sequences(candidate / "AA.fsa", assignments)
    owners = {entry.id: species for species, rows in sequences.items() for entry in rows}
    mapped, reference, ambiguous = map_cohort(cohort["pairs"], assignments, owners, status)
    if len(mapped) != 240:
        raise ValueError("Changed full experimental cohort")
    output.mkdir(parents=True, exist_ok=False)
    inputs = output / "input"
    inputs.mkdir()
    for species, rows in sequences.items():
        SeqIO.write(sorted(rows, key=lambda row: row.id), inputs / (species + ".fasta"), "fasta")
    reference_path = output / "reference_groups.json"
    with reference_path.open("x") as handle:
        json.dump(reference, handle, indent=2, sort_keys=True)
        handle.write("\n")
    report = {"status": "wgd_complete_proteomes_and_reference_mapping_prepared", "prediction_outcomes_read": False,
              "inference_authorized": False, "source": record(__file__), "source_snapshot": sources,
              "helper_sources": [record(Path(__file__).with_name(name)) for name in
                                 ("audit_ygob_overlap.py", "prepare_ygob_validation.py")],
              "cohort": record(cohort_path), "protocol": record(protocol), "species": list(SPECIES),
              "proteins_by_species": {species: len(rows) for species, rows in sequences.items()},
              "proteins": len(owners), "inputs": [record(p) for p in sorted(inputs.glob("*.fasta"))],
              "reference": record(reference_path), "reference_groups": len(reference),
              "reference_genes": sum(map(len, reference.values())), "ambiguous_source_pillars": ambiguous,
              "input_exclusions": excluded, "cohort_pairs": mapped,
              "eligible_counts": {"all": len(mapped), "split": sum(r["split_eligible"] for r in mapped),
                                  "homolog_support": sum(r["reference_eligible"] for r in mapped)},
              "prespecified_examples": cohort["prespecified_examples"],
              "limitations": ["Complete available ON protein sets for the four species, not cohort-only inputs.",
                              "Species assignments use documented YGOB columns; pillar labels are not passed to inference.",
                              "One terminal stop is removed; internal-stop proteins excluded and reported; other accepted residues preserved.",
                              "All240 experimental pairs retained; mapping eligibility is not a method outcome.",
                              "Pillars are homology groups, not copy-specific A/B orthology; ambiguous pillars excluded from reference only.",
                              "Chronology/source audit and contrast freeze remain required before inference."]}
    with (output / "manifest.json").open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("candidate", "cohort", "audit", "protocol", "output", "summary"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.summary.exists():
        raise FileExistsError(args.summary)
    result = prepare(args.candidate.resolve(), args.cohort.resolve(), args.audit.resolve(),
                     args.protocol.resolve(), args.output.resolve())
    with args.summary.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")

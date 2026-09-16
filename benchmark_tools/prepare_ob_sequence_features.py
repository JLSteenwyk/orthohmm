"""Freeze prediction-independent protein features for OrthoBench error analysis."""

import argparse
from collections import Counter
import csv
import json
import math
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH
from benchmark_tools.run_simulation_methods import read_frozen

CANONICAL = "ACDEFGHIKLMNPQRSTVWY"
LETTERS = frozenset("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
FIELDS = ("gene", "proteome", "raw_length", "residue_length", "canonical_length", "noncanonical_letter_count",
          "canonical_fraction", "normalized_entropy20", "largest_canonical_fraction", "stop_symbols",
          "gap_symbols", "other_symbols", "short_lt100", "composition_concentrated")


def sequence_features(sequence):
    if not sequence.isascii():
        raise ValueError("Non-ASCII protein sequence")
    counts = Counter(sequence.upper())
    canonical = sum(counts[letter] for letter in CANONICAL)
    residues = sum(counts[letter] for letter in LETTERS)
    fraction = canonical / residues if residues else None
    entropy = -sum((counts[letter] / canonical) * math.log2(counts[letter] / canonical)
                   for letter in CANONICAL if counts[letter]) / math.log2(20) if canonical else None
    concentrated = (entropy < .8) if canonical >= 20 and fraction >= .9 else None
    return {"raw_length": len(sequence), "residue_length": residues, "canonical_length": canonical,
            "noncanonical_letter_count": residues - canonical, "canonical_fraction": fraction,
            "normalized_entropy20": entropy,
            "largest_canonical_fraction": max((counts[letter] for letter in CANONICAL), default=0) / canonical if canonical else None,
            "stop_symbols": counts["*"], "gap_symbols": counts["-"] + counts["."],
            "other_symbols": sum(count for letter, count in counts.items() if letter not in LETTERS | {"*", "-", "."}),
            "short_lt100": residues < 100 if residues else None,
            "composition_concentrated": concentrated}


def prepare(root, output):
    root, output = root.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(output)
    source = root / "benchmark_tools/results/orthobench_factorial_prepared_20260916.json"
    prepared = read_frozen(source, PREPARED_HASH)
    inputs = prepared["fasta_inputs"]
    if len(inputs) != 12 or len({item["path"] for item in inputs}) != 12:
        raise ValueError("Expected twelve unique frozen proteome files")
    for item in inputs:
        verify_file(Path(item["path"]), item)
    output.mkdir(parents=True)
    table = output / "sequence_features.tsv"
    seen, counts, proteomes = set(), Counter({flag + "_" + state: 0
        for flag in ("short_lt100", "composition_concentrated") for state in ("true", "false", "missing")}), {}
    with table.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        for item in sorted(inputs, key=lambda record: record["path"]):
            path = Path(item["path"])
            if path.name in proteomes:
                raise ValueError("Ambiguous proteome basename")
            proteomes[path.name] = 0
            for protein in SeqIO.parse(path, "fasta"):
                if not protein.id or protein.id in seen:
                    raise ValueError("Empty or duplicate protein identifier")
                seen.add(protein.id)
                features = sequence_features(str(protein.seq))
                writer.writerow({"gene": protein.id, "proteome": path.name,
                                 **{key: "NA" if value is None else value for key, value in features.items()}})
                proteomes[path.name] += 1
                for key in ("short_lt100", "composition_concentrated"):
                    counts[key + ("_missing" if features[key] is None else "_true" if features[key] else "_false")] += 1
                for key in ("stop_symbols", "gap_symbols", "other_symbols", "noncanonical_letter_count"):
                    counts["proteins_with_" + key] += int(features[key] > 0)
                counts["empty_residue_sequences"] += int(features["residue_length"] == 0)
    if len(seen) != 251378 or any(count == 0 for count in proteomes.values()):
        raise ValueError("Incomplete frozen protein universe")
    for item in inputs:
        verify_file(Path(item["path"]), item)
    report = {"status": "sequence_features_prepared_unscored", "accuracy_evaluated": False,
              "job_id": os.environ.get("SLURM_JOB_ID"), "python": sys.version,
              "prediction_files_read": [], "reference_labels_read": [], "genes": len(seen), "proteomes": proteomes,
              "inputs": inputs, "input_manifest": file_provenance(source), "source": file_provenance(Path(__file__)),
              "table": file_provenance(table), "feature_counts": dict(counts), "missing_value": "NA",
              "definitions": {
                  "residue_length": "Count of ASCII letters A-Z after uppercase conversion; stops, gaps and other symbols excluded.",
                  "canonical_length": "Count of ACDEFGHIKLMNPQRSTVWY; noncanonical letters are separately counted, not reassigned.",
                  "normalized_entropy20": "Shannon entropy of canonical-residue frequencies divided by log2(20); NA when no canonical residues.",
                  "short_lt100": "Residue length <100, NA at zero residues; descriptive short-sequence flag, not a fragmentation diagnosis.",
                  "composition_concentrated": "Normalized entropy <0.8, evaluated only with >=20 canonical residues and canonical fraction >=0.9; otherwise NA."},
              "limitations": [
                  "Thresholds are exploratory descriptors frozen before joining error outcomes; not validated biological classifiers.",
                  "Global entropy does not identify local low-complexity tracts or domain architecture.",
                  "Sequence length alone does not establish fragmentation or completeness.",
                  "These features do not supply divergence, duplication history, domain annotation or biological ground truth.",
                  "This is a development-exposed input inventory, not independent validation or a scored error analysis."]}
    (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.root, args.output)


if __name__ == "__main__":
    main()

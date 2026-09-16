"""Validate and convert native root HOGs without inspecting reference labels.

This partition gate supplements, but does not replace, execution provenance
and full native reconciliation validation.
"""

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.score_ygob_groups import read_predictions, membership
from benchmark_tools.analyze_phylogeny_changes import split_summary
from benchmark_tools.orthobench_stage_diagnostics import file_provenance


def validate_partition(candidate_path, root_path, expected_genes, counts):
    with candidate_path.open() as handle:
        candidates = [line.split() for line in handle if line.strip()]
    candidate_groups = {f"Family{i:07d}": genes for i, genes in enumerate(candidates)}
    candidate_index = membership(candidate_groups)
    if set(candidate_index) != expected_genes:
        raise ValueError("Candidate partition differs from complete FASTA gene universe")
    predictions = read_predictions(root_path, "root_hogs")
    root_index = membership(predictions)
    if set(root_index) != expected_genes:
        raise ValueError("Root partition loses or adds input genes")
    if list(predictions) != [f"RootHOG{i:07d}" for i in range(len(predictions))]:
        raise ValueError("Root-HOG identifiers differ from native sequential order")
    with root_path.open() as handle:
        sources = [row["source_family"] for row in csv.DictReader(handle, delimiter="\t")]
    if any(source not in candidate_groups for source in sources):
        raise ValueError("Unknown or noncanonical source family")
    groups = [set(genes) for genes in predictions.values()]
    summary = split_summary([set(genes) for genes in candidates], groups, sources, set())
    for key in ("candidate_families", "root_hogs"):
        if counts.get(key) != summary[key]:
            raise ValueError(f"Native summary count mismatch: {key}")
    if counts.get("bypassed_families", -1) + counts.get("reconciled_families", -1) != len(candidates):
        raise ValueError("Native family completion counts do not cover candidates")
    summary.pop("reference_bearing_split_source_families")
    return groups, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--root-hogs", type=Path, required=True)
    parser.add_argument("--fasta-directory", type=Path, required=True)
    parser.add_argument("--native-summary", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--groups-output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if args.groups_output.exists() or args.report.exists():
        raise FileExistsError("Conversion outputs already exist")
    fastas = sorted(p for p in args.fasta_directory.iterdir() if p.suffix.lower() in {".fa", ".faa", ".fasta", ".fsa"})
    files = [args.candidate, args.root_hogs, args.native_summary, args.metrics, *fastas]
    before = [file_provenance(p) for p in files]
    genes = set()
    for path in fastas:
        for record in SeqIO.parse(path, "fasta"):
            if record.id in genes:
                raise ValueError("Duplicate input FASTA gene identifier")
            genes.add(record.id)
    if not genes:
        raise ValueError("Empty FASTA input universe")
    counts = json.loads(args.native_summary.read_text())
    metrics = json.loads(args.metrics.read_text())
    if metrics.get("status") != "complete" or metrics.get("counts") != {k: v for k, v in counts.items() if k != "schema_version"}:
        raise ValueError("Replay metrics disagree with native completion summary")
    groups, summary = validate_partition(args.candidate, args.root_hogs, genes, counts)
    if [file_provenance(p) for p in files] != before:
        raise ValueError("Inputs changed during partition validation")
    with args.groups_output.open("x") as handle:
        for group in groups:
            handle.write(" ".join(sorted(group)) + "\n")
    report = {"schema_version": 1, "status": "partition_conversion_verified",
              "input_files": before, "converter": file_provenance(Path(__file__)),
              "summary": summary, "converted_groups": file_provenance(args.groups_output),
              "accuracy_evaluated": False, "scoring_admitted": False,
              "limitations": "Execution provenance and full native reconciliation gates remain separate; no reference labels read."}
    with args.report.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()

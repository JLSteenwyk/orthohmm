"""Assemble the frozen YGOB evaluation only after fresh label-blind admission."""

import argparse
from itertools import combinations
import json
from pathlib import Path
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.report_ygob_validation import assemble_report, markdown, read_checkpoint
from benchmark_tools.score_ygob_groups import read_predictions
from benchmark_tools.verify_ygob_native_outputs import verify as verify_native
from benchmark_tools.verify_ygob_overlap_screen import verify as verify_overlap


def enumerated_counts(predictions, references):
    owners = {gene: group for group, genes in references.items() for gene in genes}
    tp = fp = 0
    for genes in predictions.values():
        for a, b in combinations([g for g in genes if g in owners], 2):
            if owners[a] == owners[b]:
                tp += 1
            else:
                fp += 1
    true_pairs = sum(1 for genes in references.values() for _ in combinations(genes, 2))
    return {"tp": tp, "fp": fp, "fn": true_pairs - tp}


def write(path, report):
    with path.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    native, overlap = verify_native(root), verify_overlap(root)
    output.mkdir(parents=True)
    write(output / "admission.json", {"native": native, "overlap": overlap, "accuracy_evaluated": False})
    references = json.loads(Path(native["reference_reconstruction"]["reference"]["path"]).read_text())
    universe = {r.id for item in native["reference_reconstruction"]["prepared_inputs"]
                for r in SeqIO.parse(item["path"], "fasta")}
    predictions = {}
    for method, evidence in native["native_groups"].items():
        path = Path(evidence["native"]["path"])
        verify_file(path, evidence["native"])
        if evidence["format"] == "native_mcl_checkpoint":
            ids = Path(evidence["sequence_ids"]["path"])
            verify_file(ids, evidence["sequence_ids"])
            predictions[method] = read_checkpoint(path, ids, universe)
        else:
            predictions[method] = read_predictions(path, evidence["format"])
    result = assemble_report(predictions, references, universe)
    enumerated = {method: enumerated_counts(groups, references) for method, groups in predictions.items()}
    for method, counts in enumerated.items():
        if counts != result["scores"][method]["counts"]:
            raise ValueError("Independent pair enumeration disagrees with scorer: " + method)
    for evidence in native["native_groups"].values():
        verify_file(Path(evidence["native"]["path"]), evidence["native"])
    reference = native["reference_reconstruction"]["reference"]
    verify_file(Path(reference["path"]), reference)
    result.update(admission=file_provenance(output / "admission.json"),
                  independent_enumerated_counts=enumerated, source=file_provenance(Path(__file__)),
                  frozen_evaluation_gates_verified=True)
    write(output / "results.json", result)
    with (output / "results.md").open("x") as handle:
        handle.write(markdown(result))
    print(output / "results.md")


if __name__ == "__main__":
    main()

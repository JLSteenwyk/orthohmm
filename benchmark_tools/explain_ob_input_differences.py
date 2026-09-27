"""Describe retained comparator input changes without repairing or promoting them."""

import argparse
import hashlib
import json
from pathlib import Path

from benchmark_tools.audit_ob_orthofinder_provenance import fasta
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def differences(original, staged, references):
    rows = []
    for gene in sorted(original.keys() & staged.keys()):
        a, b = original[gene], staged[gene]
        if a == b:
            continue
        rows.append(dict(gene=gene, original_length=len(a), staged_length=len(b),
                         original_sequence_sha256=hashlib.sha256(a.encode("ascii")).hexdigest(),
                         staged_sequence_sha256=hashlib.sha256(b.encode("ascii")).hexdigest(),
                         deletion_of_all_asterisks_explains=a.replace("*", "") == b,
                         original_asterisks=a.count("*"), staged_asterisks=b.count("*"),
                         terminal_asterisks_only=a.rstrip("*") == b,
                         reference_families=sorted(name for name, genes in references.items() if gene in genes)))
    return dict(only_original=sorted(original.keys() - staged.keys()), only_staged=sorted(staged.keys() - original.keys()),
                changed=rows, changed_count=len(rows), direct_reference_genes=sum(bool(r["reference_families"]) for r in rows))


def run(report_path, reference_path, output):
    if output.exists():
        raise FileExistsError(output)
    report_record, reference_record = record(report_path), record(reference_path)
    if reference_record["sha256"] != "660ead29c5b6ac0b8278cd1e62cdcdb0a513db81dda317e634b805e661d70ba9":
        raise ValueError("Unexpected reference inventory")
    report, reference = json.loads(report_path.read_text()), json.loads(reference_path.read_text())
    if report["status"] != "retained_ob_matrix_provenance_readback":
        raise ValueError("Unexpected input audit")
    checked = [report_record, reference_record, *report["checked_records"], *reference["inputs"]["references"]]
    for item in checked:
        check(item)
    refs = {Path(r["path"]).name: set(Path(r["path"]).read_text().splitlines()) for r in reference["inputs"]["references"]}
    if len(refs) != 70:
        raise ValueError("Incomplete reference family inventory")
    rows = []
    for method in report["rows"]:
        for item in method["staged_inputs"]:
            if item["bytes_match"]:
                continue
            diff = differences(fasta(Path(item["source"]["path"])), fasta(Path(item["staged"]["path"])), refs)
            rows.append(dict(method=method["key"], species=Path(item["staged"]["path"]).name,
                             source=item["source"], staged=item["staged"], **diff))
    for item in checked:
        check(item)
    result = dict(status="retained_input_difference_description", rows=rows, input_audit=report_record,
                  reference=reference_record, checked_records=checked, source=record(__file__),
                  changed_proteins=sum(r["changed_count"] for r in rows),
                  direct_reference_genes=sum(r["direct_reference_genes"] for r in rows),
                  input_identity_restored=False, prediction_scores_changed=False, publication_ready=False,
                  limitations=["Reference exposure uses the full 70-family memberships, before low-certainty exclusions.",
                               "A sequence transformation consistent with observed changes does not identify who performed it or when.",
                               "No counterfactual inference or accuracy-effect bound; changes outside references can alter false positives and clustering."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("report", "reference", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.report.resolve(), args.reference.resolve(), args.output.absolute())

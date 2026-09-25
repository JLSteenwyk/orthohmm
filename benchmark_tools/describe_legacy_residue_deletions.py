"""Describe reviewed native O deletions without admitting a database or search."""

import argparse
from itertools import zip_longest
import json
from pathlib import Path
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_orthomcl_database import compare_dump
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def describe(fasta, dump, reviewed_content):
    content = compare_dump(fasta, dump)
    if content != reviewed_content or not content["differences"]:
        raise ValueError("Database comparison differs from reviewed evidence or has no differences")
    expected = {row["id"]: row for row in content["differences"]}
    changes = []
    with fasta.open() as source, dump.open() as native:
        for ordinal, (a, b) in enumerate(zip_longest(
                SeqIO.parse(source, "fasta"), SeqIO.parse(native, "fasta"))):
            if a is None or b is None:
                raise ValueError("Sequence count changed during description")
            if a.id not in expected:
                continue
            before, after = str(a.seq), str(b.seq)
            positions = [i + 1 for i, residue in enumerate(before) if residue == "O"]
            if (not positions or before.replace("O", "") != after
                    or ordinal != expected[a.id]["ordinal"]):
                raise ValueError("Difference is not solely the reviewed O deletion")
            changes.append({**expected[a.id], "deleted_residue": "O",
                            "positions_one_based": positions})
    if len(changes) != len(expected):
        raise ValueError("Reviewed differences not exhausted")
    return dict(status="reviewed_native_O_deletions_described_not_admitted",
                content=content, transformations=changes,
                exact_sequence_parity=False, database_admitted=False,
                search_admitted=False, publication_ready=False,
                limitations=["Describes the extracted database, not query-side parsing.",
                             "Does not establish reference exposure or counterfactual score impact.",
                             "Does not authorize downstream inference or waive exact parity."])


def audit(report_path, report_sha256, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    report_record = record(report_path)
    if report_record["sha256"] != report_sha256:
        raise ValueError("Reviewed database report digest differs")
    report = json.loads(report_path.read_text())
    if report["status"] != "database_sequence_differences_require_review":
        raise ValueError("Expected a database difference report")
    fasta = Path(report["command"][report["command"].index("-d") + 1])
    dump = report_path.parent / "database.fasta"
    checked = [report_record, *report["checked_records"], *report["outputs"],
               record(__file__), record(Path(__file__).with_name("audit_orthomcl_database.py"))]
    if not any(row["path"] == str(fasta) for row in checked) or not any(
            row["path"] == str(dump) for row in checked):
        raise ValueError("Missing source or dump identity")
    for item in checked:
        check(item)
    result = describe(fasta, dump, report["content"])
    for item in checked:
        check(item)
    result["checked_records"] = checked
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database-report", type=Path, required=True)
    parser.add_argument("--database-report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.database_report.resolve(), args.database_report_sha256, args.output.absolute())

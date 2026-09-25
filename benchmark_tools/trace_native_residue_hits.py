"""Extract search evidence for the seven reviewed native residue deletions."""

import argparse
import hashlib
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.reviewed_legacy_database import REVIEW_SHA, validate_content
from benchmark_tools.run_simulation_methods import read_frozen


def scan(table, targets, destination, expected_sha, expected_rows):
    if not targets or any(not isinstance(t, str) or not t for t in targets):
        raise ValueError("Require nonempty target identifiers")
    if type(expected_rows) is not int or expected_rows < 1:
        raise ValueError("Require positive expected row count")
    encoded = {t.encode(): t for t in targets}
    counts = {t: dict(outgoing_hsps=0, incoming_hsps=0, self_hsps=0) for t in targets}
    outgoing, incoming = {t: set() for t in targets}, {t: set() for t in targets}
    digest = hashlib.sha256()
    rows = retained = size = 0
    # Hash the same bytes that are inspected, avoiding a separate full-table pass.
    with table.open("rb") as source, destination.open("xb") as selected:
        for line in source:
            digest.update(line)
            size += len(line)
            rows += 1
            fields = line.split(b"\t", 2)
            if len(fields) != 3:
                raise ValueError(f"Invalid BLAST identifier columns at row {rows}")
            query, subject = fields[:2]
            if query not in encoded and subject not in encoded:
                continue
            if len(line.rstrip(b"\r\n").split(b"\t")) != 12:
                raise ValueError(f"Invalid selected BLAST row {rows}")
            selected.write(line)
            retained += 1
            if query in encoded:
                key = encoded[query]
                counts[key]["outgoing_hsps"] += 1
                outgoing[key].add(subject.decode())
                counts[key]["self_hsps"] += query == subject
            if subject in encoded:
                key = encoded[subject]
                counts[key]["incoming_hsps"] += 1
                incoming[key].add(query.decode())
    if digest.hexdigest() != expected_sha or rows != expected_rows:
        raise ValueError("BLAST digest or row count differs; retained subset is unverified")
    return dict(table_sha256=digest.hexdigest(), table_bytes=size, table_rows=rows,
        retained_hsp_rows=retained, proteins=[dict(gene=t, **counts[t],
            outgoing_partners=sorted(outgoing[t]), incoming_partners=sorted(incoming[t]),
            reciprocal_partners=sorted(outgoing[t] & incoming[t])) for t in sorted(targets)])


def run(table, review_path, expected_sha, expected_rows, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    review = read_frozen(review_path, REVIEW_SHA)
    validate_content(dict(status="database_sequence_differences_require_review",
                          content=review["content"]), review)
    targets = {item["id"] for item in review["transformations"]}
    if len(targets) != 7:
        raise ValueError("Require all seven reviewed proteins")
    output.mkdir(parents=True)
    result = dict(status="selected_hit_scan_running", accuracy_admitted=False,
        publication_ready=False, source=record(Path(__file__)), review=record(review_path),
        table_path=str(table), expected_sha256=expected_sha, expected_rows=expected_rows,
        limitations=["Descriptive search trace, not independent search admission.",
                     "HSP rows are not distinct pairs or final ortholog predictions.",
                     "Self hits are included in partner lists and reciprocal counts.",
                     "No final-group membership or counterfactual accuracy effect is established."])
    report = output / "report.json"
    try:
        result["content"] = scan(table, targets, output / "selected.m8", expected_sha, expected_rows)
        result["subset"] = record(output / "selected.m8")
        result["status"] = "selected_hit_trace_verified_against_supplied_table_identity"
    except BaseException as error:
        result.update(status="selected_hit_scan_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("table", "review", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--rows", type=int, required=True)
    args = parser.parse_args()
    run(args.table.resolve(), args.review.resolve(), args.sha256, args.rows, args.output.absolute())

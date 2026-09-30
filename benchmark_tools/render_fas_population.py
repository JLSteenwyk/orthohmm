"""Display complete FAS recounts; never relabel conditional bounds as intervals."""

import argparse
import json
import math
from pathlib import Path

from benchmark_tools.map_corrected_vgnc_blocks import MANIFEST, MANIFEST_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def validate(report, manifest):
    methods, rows = manifest["methods"], report["methods"]
    if (report.get("status") != "retained_fas_eligible_populations_recounted"
            or report.get("uncertainty_admitted") is not False
            or report.get("benchmark_scores_changed") is not False
            or report.get("publication_ready") is not False
            or len(rows) != 8 or len(methods) != 8
            or len({r["method"] for r in rows}) != 8
            or [r["method"] for r in rows] != [m["key"] for m in methods]):
        raise ValueError("Require the complete unchanged eight-method recount")
    for row, method in zip(rows, methods):
        if (row["native_counts_match"] is not True
                or row["saved_lookup_strata_and_values_match"] is not True
                or row["database_historically_hash_bound"] is not False):
            raise ValueError("Invalid recount/sample or historical database claim")
        fields = ("distinct_query_pairs", "skipped_alias_pairs", "precomputed", "missing", "unannotated", "eligible_pairs")
        if any(type(row[k]) is not int or row[k] < 0 for k in fields):
            raise ValueError("Population counts must be nonnegative integers")
        p, m, n, s = row["precomputed"], row["missing"], row["eligible_pairs"], row["precomputed_score_sum"]
        if (n <= 0 or n != p + m or n != method["details"]["FAS"]["assessed_relations"]
                or row["distinct_query_pairs"] != n + row["unannotated"] + row["skipped_alias_pairs"]
                or not isinstance(s, (int, float)) or not math.isfinite(s) or not 0 <= s <= p):
            raise ValueError("Counts or score sum inconsistent")
        if any(row["native_logged_counts"][k] != row[k] for k in ("precomputed", "missing", "unannotated")):
            raise ValueError("Native count equality is not reproduced")
        expected = [s / n, (s + m) / n]
        observed = row["hypothetical_full_mean_bounds"]
        if len(observed) != 2 or any(a != b for a, b in zip(observed, expected)):
            raise ValueError("Completion-bound arithmetic differs")
    return rows


def render(report, manifest):
    rows = validate(report, manifest)
    lines = ["# Retained FAS Eligible-Population Recount", "",
        "Generated from a complete eight-method recount. Completion bounds assume",
        "uncomputed eligible scores have hypothetical values in [0,1]. They are not",
        "confidence intervals or replacement native FAS scores. Database hashes are",
        "new retained identities, not independent historical checksum bindings.", "",
        "| Method | Eligible pairs | Precomputed | Uncomputed | Precomputed mean | Hypothetical full-mean bound |",
        "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for row in rows:
        mean = row["precomputed_score_sum"] / row["precomputed"] if row["precomputed"] else None
        label = next(m["label"] for m in manifest["methods"] if m["key"] == row["method"])
        low, high = row["hypothetical_full_mean_bounds"]
        lines.append(f"| {label} | {row['eligible_pairs']:,} | {row['precomputed']:,} | "
                     f"{row['missing']:,} | {'NA' if mean is None else format(mean, '.6f')} | {low:.6f} - {high:.6f} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    report_ref, manifest_ref = record(args.report), record(args.repo.resolve() / MANIFEST)
    if report_ref["sha256"] != args.report_sha256 or manifest_ref["sha256"] != MANIFEST_SHA:
        raise ValueError("Report or corrected manifest digest differs")
    report = json.loads(Path(report_ref["path"]).read_text())
    manifest = json.loads(Path(manifest_ref["path"]).read_text())
    for pin in report["checked_records"]:
        check(pin)
    text = render(report, manifest)
    check(report_ref)
    check(manifest_ref)
    with args.output.open("x") as stream:
        stream.write(text)

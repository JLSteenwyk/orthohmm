"""Check annotation-path panel completeness without claiming independent graph validation."""

import argparse
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def validate(rows, protein_ids, options):
    ids = [r["protein"] for r in rows]
    if len(set(ids)) != len(ids) or set(ids) != set(protein_ids):
        raise ValueError("Missing, extra or duplicate protein records")
    limit = options["paths_limit"]
    if type(limit) is not int or limit != 10 ** 15:
        raise ValueError("Unexpected path limit")
    for row in rows:
        count = row["paths"]
        instances = row["linearized_instances"]
        if type(count) is not int or count < 0 or type(instances) is not int or instances < 0:
            raise ValueError("Invalid path or feature count")
        if type(row["rejected_by_path_limit"]) is not bool or row["rejected_by_path_limit"] != (count > limit):
            raise ValueError("Exclusion flag differs from strict native limit")
    return dict(proteins=len(rows), rejected=sum(r["rejected_by_path_limit"] for r in rows),
                maximum_paths=max((r["paths"] for r in rows), default=0))


def review(path):
    source = record(path)
    panel = json.loads(path.read_text())
    if panel["status"] != "retained_annotation_path_panel_complete":
        raise ValueError("Require complete panel")
    files = panel["files"]
    expected = [r for r in panel["checked_inputs"] if "/fas_annotations/" in r["path"]]
    if len(files) != len(expected) or sorted(r["annotation"]["path"] for r in files) != sorted(r["path"] for r in expected):
        raise ValueError("Panel differs from frozen annotation inventory")
    all_ids, repeated, rejected = set(), set(), set()
    summaries = []
    checked = [source, *panel["checked_inputs"]]
    for ref in checked:
        check(ref)
    for item in files:
        for ref in (item["annotation"], item["result"]):
            check(ref)
            checked.append(ref)
        data = json.loads(Path(item["result"]["path"]).read_text())
        annotation = json.loads(Path(item["annotation"]["path"]).read_text())
        if data["annotation"] != item["annotation"] or data["options"] != item["options"]:
            raise ValueError("Per-file provenance differs")
        summary = validate(data["proteins"], annotation["feature"], data["options"])
        if summary != item["summary"] or summary != data["summary"]:
            raise ValueError("Per-file summary differs from rows")
        ids = set(annotation["feature"])
        repeated.update(all_ids & ids)
        all_ids.update(ids)
        rejected.update(r["protein"] for r in data["proteins"] if r["rejected_by_path_limit"])
        summaries.append(summary)
    totals = dict(annotation_files=len(files), protein_records=sum(s["proteins"] for s in summaries),
                  rejected_records=sum(s["rejected"] for s in summaries))
    if totals != panel["summary"]:
        raise ValueError("Panel summary differs from rows")
    for ref in checked:
        check(ref)
    return dict(status="native_path_panel_record_completeness_verified", source=record(__file__),
        panel=source, summary=totals, unique_protein_identifiers=len(all_ids),
        identifiers_in_multiple_annotation_files=len(repeated), unique_rejected_identifiers=len(rejected),
        historical_omissions_attributed=False,
        limitations=["Verifies complete identifier coverage and summary arithmetic, not an independent graph algorithm.",
            "Counts are not historical sampled-pair omissions or evidence of comparison bias."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = review(args.panel)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")

"""Prospective NA serialization correction; original failed reader stays frozen."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

from benchmark_tools import readback_controlled_fragment_stages as previous


PREVIOUS_SHA = "f353ef1d9603a1bf1f4ee50db015d0f23fb24fd07e423590c74590e00a643952"
FIELDS = ("case_id", "method", "seed", "fragment_endpoints", "category", "truth", "arm", "predicted",
          "search_status", "hit_forward", "hit_reverse", "graph_direct", "graph_connected", "same_seed_group",
          "same_candidate", "same_root_hog", "pair_event", "membership_filter_active", "observed_location")


def serialized(row):
    return {k: "NA" if row[k] is None else str(row[k]) for k in FIELDS}


def validate_table(path, expected):
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        previous.need(tuple(reader.fieldnames or ()) == FIELDS, "Changed stage table header")
        rows = list(reader)
    previous.need(rows == [serialized(r) for r in expected], "Stage TSV value mismatch")
    return len(rows)


def verify(report_path, expected_sha):
    source = previous.pin(previous.__file__)
    previous.need(source["sha256"] == PREVIOUS_SHA, "Changed original independent reader")
    ref = previous.pin(report_path)
    previous.need(ref["sha256"] == expected_sha, "Changed stage report")
    report = json.loads(Path(report_path).read_text())
    previous.need(report["schema"] == "controlled_fragment_stage_trace_v1" and report["status"] == "retained_stages_verified"
                  and report["new_inference_or_scoring"] is False and report["new_bootstrap_draws"] == 0
                  and all(report[k] is False for k in ("publication_ready", "uncertainty_admitted", "scientific_timings_admitted", "independent_confirmation")), "Changed stage scope")
    selected = json.loads(previous.check(report["selection"]).read_text())
    checked = [previous.pin(previous.check(r)) for r in report["checked_inputs"]]
    bindings = {b["seed"]: b["arms"] for b in selected["bindings"]}
    query_sets = {}
    for case in selected["cases"]:
        for arm in previous.ARMS:
            query_sets.setdefault((case["method"], case["seed"], arm), set()).add((case["gene_a"], case["gene_b"]))
    contexts = {key: previous.read_context(bindings[key[1]][key[2]][key[0]], key[0], queries) for key, queries in query_sets.items()}
    expected, counts = [], Counter()
    previous.need(len(report["cases"]) == len(selected["cases"]), "Changed case count")
    for case, original in zip(report["cases"], selected["cases"]):
        previous.need({k:v for k,v in case.items() if k != "stages"} == original, "Changed selected identity")
        for arm in previous.ARMS:
            row = previous.observation(case["method"], contexts[case["method"], case["seed"], arm], (case["gene_a"], case["gene_b"]))
            previous.need(row == case["stages"][arm], "Native stage mismatch: " + case["case_id"] + "/" + arm)
            previous.need(row["predicted"] is case["comparator_predictions"][arm][case["method"]], "Comparator mismatch")
            expected.append(dict(case_id=case["case_id"], method=case["method"], seed=case["seed"], fragment_endpoints=case["fragment_endpoints"],
                                 category=case["category"], truth=case["truth"], arm=arm, **row))
            counts[case["method"], case["category"], row["observed_location"], arm] += 1
    summary = [dict(method=m, category=c, observed_location=l, arm=a, representatives=n) for (m,c,l,a),n in sorted(counts.items())]
    previous.need(summary == report["summary"], "Changed stage summary")
    rows = validate_table(Path(report_path).with_name("stages.tsv"), expected)
    previous.need(rows == report["stage_rows"], "Changed stage row count")
    for r in checked:
        previous.check(r)
    return dict(status="independent_native_stages_and_na_table_verified", report=ref, stage_rows=rows,
                selected_cases=len(report["cases"]), native_contexts=len(contexts), checked_input_files=len(checked),
                summary=summary, source=previous.pin(__file__), previous_reader=source,
                compatibility_mapping={"JSON null": "TSV NA", "JSON false": "TSV False", "JSON true": "TSV True"},
                original_report_modified=False, new_inference_or_scoring=False, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    previous.need(not args.output.exists() and not args.output.is_symlink(), "Existing readback output")
    result = verify(args.report, args.sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k:result[k] for k in ("status", "stage_rows", "selected_cases", "native_contexts")}))

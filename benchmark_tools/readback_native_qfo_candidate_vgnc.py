"""Independently read candidate VGNC rows and cached baseline pair states."""

import argparse
from collections import Counter
import csv
from fractions import Fraction
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import readback_native_qfo_vgnc_blocks as independent

ROOT = Path(__file__).resolve().parent.parent
CELLS = ("p0_c0_r0", "p0_c1_r0")


def state(value):
    tokens = set() if value == "not_scored" else set(value.split("+"))
    canonical = "+".join(c for c in independent.CATEGORIES if c in tokens) or "not_scored"
    independent.need(tokens <= set(independent.CATEGORIES) and canonical == value, "Invalid baseline state")
    return tokens


def metrics(counts):
    tp, fp, fn = (counts[c] for c in independent.CATEGORIES)
    independent.need(tp + fp > 0 and tp + fn > 0, "Undefined pooled ratio")
    return dict(precision=Fraction(tp, tp + fp), recall=Fraction(tp, tp + fn),
                f1=Fraction(2 * tp, 2 * tp + fp + fn))


def review(path, digest):
    ref = independent.identity(path)
    independent.need(ref["sha256"] == digest, "Changed candidate decomposition report")
    report = independent.read_json(ref)
    independent.need(report.get("schema") == "native_qfo_candidate_vgnc_decomposition_v1"
        and report.get("status") == "new_candidate_scored_rows_decomposed"
        and all(report.get(k) is False for k in ("uncertainty_admitted", "new_scoring_or_admission",
            "publication_ready", "baseline_raw_audit_repeated", "prediction_edges_requeried"))
        and report.get("failed_r1_timing_remains_ineligible") is True
        and report.get("full_candidate_database_hash_checked") is True, "Changed decomposition scope")
    checked = report["checked_records"]
    for record in checked:
        independent.verify(record)
    primary = independent.identity(Path(__file__).with_name("export_native_qfo_candidate_vgnc.py"))
    independent.need(report["source"] == primary and independent.identity(independent.__file__) in checked,
                     "Primary or independent-kernel source differs")
    for key in ("source", "protocol", "snapshot", "baseline_report", "baseline_readback",
                "reference_record", "baseline_pair_table"):
        independent.need(report[key] in checked, "Unbound input: " + key)
    snapshot, previous, readback = [independent.read_json(report[k]) for k in
        ("snapshot", "baseline_report", "baseline_readback")]
    independent.need(readback["report"] == report["baseline_report"]
        and readback["schema"] == "native_qfo_vgnc_blocks_readback_v1"
        and readback["status"] == "complete_native_decomposition_independently_verified"
        and all(readback.get(k) is False for k in ("uncertainty_admitted", "new_scoring_or_admission", "publication_ready"))
        and readback.get("failed_r1_timing_remains_ineligible") is True
        and readback["source"] == independent.identity(independent.__file__)
        and previous["schema"] == "native_qfo_vgnc_blocks_v1"
        and previous["status"] == "native_scored_rows_decomposed"
        and previous["source"] == independent.identity(Path(__file__).with_name("export_native_qfo_vgnc_blocks.py"))
        and all(previous.get(k) is False for k in ("uncertainty_admitted", "new_scoring_or_admission", "publication_ready"))
        and previous.get("failed_r1_timing_remains_ineligible") is True
        and previous["methods"][0] == report["baseline"]
        and previous["transition_table"] == report["baseline_pair_table"], "Cached baseline linkage differs")
    independent.need(snapshot["schema"] == "native_qfo_scientific_reporting_snapshot_v1"
        and snapshot["source"] in checked and snapshot["plan"] in checked
        and all(snapshot.get(k) is False for k in ("new_scoring_or_admission", "publication_ready",
            "recovered_inference_resources_admitted")), "Snapshot scope/source differs")
    rows = [r for r in snapshot["rows"] if r.get("accuracy_admitted") is True]
    independent.need([r["index"] for r in snapshot["rows"]] == list(range(6, 13))
        and [(r["index"], r["cell"]) for r in rows] == [(6, CELLS[0]), (7, "p0_c0_r1"), (8, CELLS[1])]
        and rows[0]["prediction_semantics"] == rows[2]["prediction_semantics"]
        and rows[1].get("timing_eligible") is False and rows[1].get("timing_admitted") is False,
        "Cohort or failed timing differs")
    for method, row in zip((report["baseline"], report["candidate"]), (rows[0], rows[2])):
        independent.need(all(method[k] == row[k] for k in ("index", "cell", "participant", "admission",
            "prediction_semantics", "measurement_status")) and row["admission"] in snapshot["evidence"],
            "Snapshot method linkage differs")
    assertions, labels, mapping, summary = independent.reconstruct(report["reference_record"]["path"])
    independent.need(summary == previous["reference"] == report["reference"]
        and len(assertions) == previous["reference_pairs"] == report["reference_pairs"], "Reference differs")
    candidate = report["candidate"]
    admission = independent.read_json(candidate["admission"])
    execution = independent.read_json(candidate["execution"])
    independent.need(candidate["admission"] in checked and candidate["execution"] in checked
        and admission["schema"] == "full_native_factorial_qfo_admission_v1"
        and admission["status"] == "full_native_factorial_qfo_assessment_admitted"
        and admission["accuracy_admitted"] is True and admission["publication_ready"] is False
        and admission["execution_report"] == candidate["execution"]
        and candidate["execution"] in admission["checked_records"]
        and (admission["native_index"], admission["cell"], admission["participant"])
            == (8, CELLS[1], candidate["participant"])
        and execution["exit_code"] == 0 and execution["native_index"] == 8
        and execution["cell"] == CELLS[1]
        and execution["status"] == "process_succeeded_pending_independent_admission", "Original admission/execution differs")
    for key in ("database", "aggregate", "raw"):
        independent.need(candidate[key] in checked and candidate[key] in admission["checked_records"]
            and candidate[key] in execution["outputs"], "Original candidate artifact inventory differs")
    selected, mapping_digest, aliases = independent.mapped_database(candidate["database"]["path"], labels)
    independent.need(mapping_digest == previous["methods"][0]["selected_reference_mapping_sha256"]
        == candidate["selected_reference_mapping_sha256"] and aliases == candidate["alias_rows"], "Mapping or aliases differ")
    annotations = {v[0]: (v[1], v[2]) for v in selected.values()}
    asserted = {tuple(sorted((selected[a][0], selected[b][0]))) for a, b in assertions}
    cells, pairs, validation, raw_rows = independent.raw_counts(candidate["raw"]["path"], asserted, annotations, mapping)
    independent.need(validation == candidate["validation"] and validation["counts"] == candidate["counts"], "Candidate counts differ")
    independent.table(candidate["table"], ["block_left", "block_right", *independent.CATEGORIES],
        [[a, b, *[v[c] for c in independent.CATEGORIES]] for (a, b), v in sorted(cells.items())])
    independent.need(candidate["nonzero_cells"] == len(cells)
        and candidate["nonzero_cross_block_cells"] == sum(a != b for a, b in cells)
        and candidate["within_block_false_positives"] == sum(v["FP"] for (a, b), v in cells.items() if a == b)
        and candidate["cross_block_false_positives"] == sum(v["FP"] for (a, b), v in cells.items() if a != b),
        "Sparse block summary differs")
    baseline_pairs, seen, counts, original_transitions = {}, set(), Counter(), Counter()
    with Path(report["baseline_pair_table"]["path"]).open(newline="") as stream:
        table = csv.DictReader(stream, delimiter="\t")
        independent.need(table.fieldnames == ["protein_left", "protein_right", "block_left", "block_right",
                                             "p0_c0_r0", "p0_c0_r1"], "Cached pair header differs")
        for row in table:
            independent.need(None not in row and None not in row.values(), "Malformed cached pair row")
            pair = (row["protein_left"], row["protein_right"])
            independent.need(pair[0] < pair[1] and pair not in seen and all(g in annotations for g in pair),
                             "Duplicate, unordered or unmapped cached pair")
            seen.add(pair)
            independent.need([row["block_left"], row["block_right"]]
                == sorted(mapping[annotations[g][0]] for g in pair), "Cached block annotation differs")
            left, right = state(row["p0_c0_r0"]), state(row["p0_c0_r1"])
            independent.need(left or right, "Unscored cached union pair")
            if left:
                baseline_pairs[pair] = left
                counts.update(left)
            original_transitions[row["p0_c0_r0"], row["p0_c0_r1"]] += 1
    independent.need(len(seen) == previous["union_scored_pairs"]
        and [dict(r0=a, r1=b, pairs=n) for (a, b), n in sorted(original_transitions.items())] == previous["transition_counts"]
        and {c: counts[c] for c in independent.CATEGORIES} == report["baseline"]["counts"], "Cached baseline counts differ")
    transitions, transition_counts = [], Counter()
    for a, b in sorted(set(baseline_pairs) | set(pairs)):
        states = ["+".join(c for c in independent.CATEGORIES if c in p.get((a, b), set())) or "not_scored"
                  for p in (baseline_pairs, pairs)]
        blocks = sorted(mapping[annotations[g][0]] for g in (a, b))
        transitions.append([a, b, *blocks, *states])
        transition_counts[tuple(states)] += 1
    independent.table(report["transition_table"], ["protein_left", "protein_right", "block_left", "block_right", *CELLS], transitions)
    independent.need(len(transitions) == report["union_scored_pairs"] and report["transition_counts"]
        == [dict(baseline=a, candidate=b, pairs=n) for (a, b), n in sorted(transition_counts.items())], "Candidate transitions differ")
    values = [metrics(m["counts"]) for m in (report["baseline"], candidate)]
    for observed, method, row in zip(values, (report["baseline"], candidate), (rows[0], rows[2])):
        for key, value in observed.items():
            independent.near(value, method["metrics"][key])
            independent.near(value, row["scores"]["VGNC"] if key == "f1" else row["endpoint_details"]["VGNC"][key])
    inline = independent.read_json(candidate["aggregate"])["datalink"]["inline_data"]
    native = inline["challenge_participants"]
    independent.need(len(native) == 1 and native[0]["participant_id"] == candidate["participant"]
        and native[0] == admission["assessment"]["endpoints"]["VGNC"]["native_participant"]
        and (inline["visualization"]["x_axis"], inline["visualization"]["y_axis"]) == ("TPR", "PPV"), "Native aggregate differs")
    independent.near(values[1]["precision"], native[0]["metric_y"])
    independent.near(values[1]["recall"], native[0]["metric_x"])
    for key in values[0]:
        independent.near(values[1][key] - values[0][key], report["differences"][key])
    for record in checked:
        independent.verify(record)
    independent.verify(ref)
    return dict(schema="native_qfo_candidate_vgnc_readback_v1", status="candidate_decomposition_independently_verified",
        source=independent.identity(__file__), independent_kernel_source=independent.identity(independent.__file__),
        report=ref, candidate_raw_rows=raw_rows, cached_baseline_union_rows=len(seen), transition_pairs=len(transitions),
        checked_input_records=len(checked), candidate_counts=validation["counts"], transition_counts=report["transition_counts"],
        primary_exporter_imported=False, baseline_raw_audit_repeated=False, uncertainty_admitted=False,
        new_scoring_or_admission=False, failed_r1_timing_remains_ineligible=True, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    independent.need(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = review(args.report, args.report_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write("\n")
    print(json.dumps(dict(status=result["status"], candidate_counts=result["candidate_counts"], transition_pairs=result["transition_pairs"])))

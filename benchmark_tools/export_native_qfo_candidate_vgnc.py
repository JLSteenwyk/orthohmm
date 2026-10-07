"""Trace newly admitted candidate VGNC rows, reusing a verified baseline."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import export_native_qfo_vgnc_blocks as old

CELLS = ("p0_c0_r0", "p0_c1_r0")
ROOT = Path(__file__).resolve().parent.parent
PROTOCOL = ROOT / "benchmark_tools/results/NATIVE_QFO_CANDIDATE_VGNC_PROTOCOL_20261006.md"
SNAPSHOT_SHA = "7916d3e23808edbb92b016b5c40ac56b9417a863d1af01a7b93a8ef5dbb53a63"
BASELINE_SHA = "23cb9fd4b1c2219239ff354c97cd7dab9758667700667f88e8424414184e4c1f"
READBACK_SHA = "6eadf57f4f3032bbc3b787c6771824f73f47baae9c643d524b7422709ee94720"


def categories(value):
    result = set() if value == "not_scored" else set(value.split("+"))
    old.require(result <= set(old.CATEGORIES) and old.status(result) == value, "Invalid cached category state")
    return result


def baseline_pairs(previous, mapping, annotations):
    old.check(previous["transition_table"])
    pairs, seen, totals, transitions = {}, set(), Counter(), Counter()
    with Path(previous["transition_table"]["path"]).open(newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        old.require(next(rows) == ["protein_left", "protein_right", "block_left", "block_right", *old.CELLS],
                    "Unexpected verified baseline table header")
        for row in rows:
            old.require(len(row) == 6, "Invalid baseline table columns")
            a, b, left, right, state, recovered = row
            pair = (a, b)
            old.require(a < b and pair not in seen and a in annotations and b in annotations,
                        "Duplicate, unordered or unmapped baseline pair")
            seen.add(pair)
            old.require([left, right] == sorted(mapping[annotations[g][0]] for g in pair),
                        "Cached baseline block annotation differs")
            values = categories(state)
            categories(recovered)
            old.require(values or recovered != "not_scored", "Unscored pair in baseline union")
            if values:
                pairs[pair] = values
                totals.update(values)
            transitions[state, recovered] += 1
    old.require(len(seen) == previous["union_scored_pairs"] and
                [dict(r0=a, r1=b, pairs=n) for (a, b), n in sorted(transitions.items())]
                == previous["transition_counts"], "Cached baseline transition summary differs")
    old.require({c: totals[c] for c in old.CATEGORIES} == previous["methods"][0]["counts"],
                "Cached baseline counts differ")
    return pairs


def export(snapshot_path, snapshot_sha, baseline_path, baseline_sha, readback_path, readback_sha, output):
    output = Path(output).absolute()
    old.require(output.resolve() == output and not output.exists() and not output.is_symlink(),
                "Require fresh direct output directory")
    checked = [old.record(__file__), old.record(PROTOCOL), old.record(old.__file__),
        old.record(ROOT / "benchmark_tools/readback_native_qfo_vgnc_blocks.py"),
        old.record(ROOT / "benchmark_tools/audit_qfo_vgnc_mapping.py"),
        old.record(ROOT / "benchmark_tools/prepare_ob_candidate_neighborhood.py")]
    refs = [old.record(path) for path in (snapshot_path, baseline_path, readback_path)]
    old.require([r["sha256"] for r in refs] == [snapshot_sha, baseline_sha, readback_sha], "Changed selected input")
    snapshot, previous, readback = [old.load(ref, checked) for ref in refs]
    old.require(snapshot.get("schema") == "native_qfo_scientific_reporting_snapshot_v1"
        and all(snapshot.get(k) is False for k in ("new_scoring_or_admission", "publication_ready",
            "recovered_inference_resources_admitted")), "Changed snapshot scope")
    old.check(snapshot["source"]); checked.append(snapshot["source"])
    old.load(snapshot["plan"], checked)
    old.require(previous.get("schema") == "native_qfo_vgnc_blocks_v1"
        and previous.get("status") == "native_scored_rows_decomposed"
        and all(previous.get(k) is False for k in ("uncertainty_admitted", "new_scoring_or_admission", "publication_ready"))
        and previous.get("failed_r1_timing_remains_ineligible") is True
        and previous["source"] == old.record(old.__file__), "Changed baseline decomposition scope/source")
    old.require(readback.get("schema") == "native_qfo_vgnc_blocks_readback_v1"
        and readback.get("status") == "complete_native_decomposition_independently_verified"
        and readback.get("report") == refs[1] and readback.get("failed_r1_timing_remains_ineligible") is True
        and all(readback.get(k) is False for k in ("uncertainty_admitted", "new_scoring_or_admission", "publication_ready"))
        and readback.get("source") == checked[3], "Baseline independent readback differs")
    rows = snapshot["rows"]
    old.require([r["index"] for r in rows] == list(range(6, 13)), "Incomplete or duplicate snapshot identities")
    admitted = [r for r in rows if r.get("accuracy_admitted") is True]
    old.require([(r["index"], r["cell"]) for r in admitted] == [(6, CELLS[0]), (7, old.CELLS[1]), (8, CELLS[1])],
                "Unexpected admitted cohort")
    baseline, recovered, candidate = admitted
    old.require(candidate["prediction_semantics"] == baseline["prediction_semantics"],
                "Candidate and baseline prediction semantics differ")
    old.require(recovered.get("timing_eligible") is False and recovered.get("timing_admitted") is False,
                "Recovered failed timing relabeled")
    original = previous["methods"][0]
    old.require([m["cell"] for m in previous["methods"]] == list(old.CELLS)
        and all(baseline[k] == original[k] for k in ("index", "cell", "participant", "prediction_semantics",
            "measurement_status", "admission")) and baseline["admission"] in snapshot["evidence"],
        "Verified baseline identity differs")
    old.close(original["metrics"], dict(precision=baseline["endpoint_details"]["VGNC"]["precision"],
        recall=baseline["endpoint_details"]["VGNC"]["recall"], f1=baseline["scores"]["VGNC"]))
    for ref in (previous["reference_record"], previous["reference_table"], previous["transition_table"]):
        old.check(ref); checked.append(ref)
    truth, labels = old.reference_data(Path(previous["reference_record"]["path"]))
    mapping, summary = old.reference_blocks(truth)
    old.require(summary == previous["reference"] and len(truth) == previous["reference_pairs"],
                "Reference reconstruction differs")
    old.require(candidate["admission"] in snapshot["evidence"], "Candidate admission unbound")
    admission = old.load(candidate["admission"], checked)
    old.require(admission.get("schema") == "full_native_factorial_qfo_admission_v1"
        and admission.get("status") == "full_native_factorial_qfo_assessment_admitted"
        and admission.get("accuracy_admitted") is True and admission.get("publication_ready") is False
        and all(admission[k] == candidate[r] for k, r in (("native_index", "index"), ("cell", "cell"),
            ("participant", "participant"))) and admission["execution_report"] in admission["checked_records"],
        "Candidate admission differs")
    execution = old.load(admission["execution_report"], checked)
    old.require(execution.get("exit_code") == 0 and execution.get("native_index") == 8
        and execution.get("cell") == CELLS[1]
        and execution.get("status") == "process_succeeded_pending_independent_admission", "Candidate execution differs")
    participant = candidate["participant"]
    artifacts = [old.unique_record(admission["checked_records"], suffix) for suffix in (
        "/other/" + participant + ".db", "/results/VGNC/VGNC.json",
        "/VGNC_" + participant.replace("_", "-").replace(" ", "-") + "_raw.txt.gz")]
    old.require(all(ref in execution["outputs"] for ref in artifacts), "Candidate artifact not originally inventoried")
    for ref in artifacts:
        old.check(ref); checked.append(ref)
    database, metric, raw = [Path(r["path"]) for r in artifacts]
    asserted, annotations, digest, aliases = old.mapped_reference(database, truth, labels)
    old.require(digest == original["selected_reference_mapping_sha256"], "Candidate reference mapping differs")
    validation = old.validate_raw(raw, asserted, annotations)
    cells, statuses = old.aggregate(raw, mapping)
    totals = {c: sum(v[c] for v in cells.values()) for c in old.CATEGORIES}
    old.require(totals == validation["counts"], "Candidate block totals differ")
    observed = old.ratios(totals)
    old.close(observed, dict(precision=candidate["endpoint_details"]["VGNC"]["precision"],
        recall=candidate["endpoint_details"]["VGNC"]["recall"], f1=candidate["scores"]["VGNC"]))
    inline = json.loads(metric.read_text())["datalink"]["inline_data"]
    native = inline["challenge_participants"]
    old.require(len(native) == 1 and native[0]["participant_id"] == participant
        and (inline["visualization"]["x_axis"], inline["visualization"]["y_axis"]) == ("TPR", "PPV")
        and native[0] == admission["assessment"]["endpoints"]["VGNC"]["native_participant"], "Native aggregate differs")
    old.close({k: observed[k] for k in ("precision", "recall")},
              dict(precision=native[0]["metric_y"], recall=native[0]["metric_x"]))
    baseline_statuses = baseline_pairs(previous, mapping, annotations)
    transition_rows, transitions = [], Counter()
    for a, b in sorted(set(baseline_statuses) | set(statuses)):
        states = (old.status(baseline_statuses.get((a, b), set())), old.status(statuses.get((a, b), set())))
        blocks = sorted(mapping[annotations[g][0]] for g in (a, b))
        transition_rows.append([a, b, *blocks, *states])
        transitions[states] += 1
    for ref in checked:
        old.check(ref)
    output.mkdir(parents=True, exist_ok=False)
    sparse, transition = output / "candidate_blocks.tsv", output / "pair_transitions.tsv"
    old.write_table(sparse, ["block_left", "block_right", *old.CATEGORIES],
        [[a, b, *[v[c] for c in old.CATEGORIES]] for (a, b), v in sorted(cells.items())])
    old.write_table(transition, ["protein_left", "protein_right", "block_left", "block_right", *CELLS], transition_rows)
    result = dict(schema="native_qfo_candidate_vgnc_decomposition_v1", status="new_candidate_scored_rows_decomposed",
        source=checked[0], protocol=checked[1], snapshot=refs[0], baseline_report=refs[1], baseline_readback=refs[2],
        reference_record=previous["reference_record"], reference=summary, reference_pairs=len(truth),
        baseline=original, baseline_pair_table=previous["transition_table"],
        candidate=dict(index=8, cell=CELLS[1], participant=participant, admission=candidate["admission"],
            execution=admission["execution_report"], database=artifacts[0], aggregate=artifacts[1], raw=artifacts[2],
            prediction_semantics=candidate["prediction_semantics"], measurement_status=candidate["measurement_status"],
            counts=totals, metrics=observed, validation=validation, alias_rows=aliases,
            selected_reference_mapping_sha256=digest, table=old.record(sparse), nonzero_cells=len(cells),
            nonzero_cross_block_cells=sum(a != b for a, b in cells),
            within_block_false_positives=sum(v["FP"] for (a, b), v in cells.items() if a == b),
            cross_block_false_positives=sum(v["FP"] for (a, b), v in cells.items() if a != b)),
        transition_table=old.record(transition), union_scored_pairs=len(transition_rows),
        transition_counts=[dict(baseline=a, candidate=b, pairs=n) for (a, b), n in sorted(transitions.items())],
        differences={k: observed[k] - original["metrics"][k] for k in observed}, checked_records=checked,
        baseline_raw_audit_repeated=False, prediction_edges_requeried=False, full_candidate_database_hash_checked=True,
        failed_r1_timing_remains_ineligible=True, uncertainty_admitted=False, new_scoring_or_admission=False,
        publication_ready=False, limitations=previous["limitations"] + [
            "Baseline uses previously independently verified complete pair states, not a new raw-data audit.",
            "Candidate aggregate outcomes were inspected before this raw decomposition was prespecified."])
    for ref in checked:
        old.check(ref)
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, default=ROOT / "benchmark_tools/results/native_qfo_scientific_scores_20261006_v2/report.json")
    parser.add_argument("--snapshot-sha256", default=SNAPSHOT_SHA)
    parser.add_argument("--baseline", type=Path, default=ROOT / "benchmark_tools/results/native_qfo_vgnc_blocks_20261006_v1/report.json")
    parser.add_argument("--baseline-sha256", default=BASELINE_SHA)
    parser.add_argument("--readback", type=Path, default=ROOT / "benchmark_tools/results/native_qfo_vgnc_blocks_readback_20261006_v1.json")
    parser.add_argument("--readback-sha256", default=READBACK_SHA)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = export(args.snapshot, args.snapshot_sha256, args.baseline, args.baseline_sha256,
                    args.readback, args.readback_sha256, args.output)
    print(json.dumps(dict(status=result["status"], counts=result["candidate"]["counts"], transitions=result["transition_counts"])))

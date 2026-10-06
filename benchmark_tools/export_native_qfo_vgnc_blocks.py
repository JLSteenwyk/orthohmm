"""Decompose admitted native VGNC rows without inventing independent units."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_vgnc_mapping import reference_data, mapped_reference, validate_raw
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

SNAPSHOT_SHA = "6b2f735ea8a44f72715e11a7f015e4328576baa74c544889cbed0cdfe70ea07b"
HISTORICAL_SHA = "01a51a62a536bcf759f2a0ef83c56196a92022b0bb8290448eacbfefa47d6005"
CELLS = ("p0_c0_r0", "p0_c0_r1")
CATEGORIES = ("TP", "FP", "FN")
PROTOCOL = Path(__file__).parent / "results/NATIVE_QFO_VGNC_BLOCK_PROTOCOL_20261006.md"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(ref, checked):
    check(ref)
    checked.append(ref)
    return json.loads(Path(ref["path"]).read_text())


def unique_record(records, suffix):
    candidates = {json.dumps(r, sort_keys=True): r for r in records if r["path"].endswith(suffix)}
    require(len(candidates) == 1, "Missing, ambiguous or conflicting inventory: " + suffix)
    return next(iter(candidates.values()))


def reference_blocks(truth):
    owners, adjacency = defaultdict(set), defaultdict(set)
    for pair, label in truth.items():
        adjacency[label]
        for protein in pair:
            owners[protein].add(label)
    for labels in owners.values():
        first = min(labels)
        for other in labels - {first}:
            adjacency[first].add(other)
            adjacency[other].add(first)
    mapping, groups = {}, []
    for label in sorted(adjacency):
        if label in mapping:
            continue
        seen, pending = {label}, [label]
        while pending:
            node = pending.pop()
            for neighbor in adjacency[node] - seen:
                seen.add(neighbor)
                pending.append(neighbor)
        group = sorted(seen)
        groups.append(group)
        mapping.update((member, group[0]) for member in group)
    return mapping, dict(family_labels=len(adjacency), reference_proteins=len(owners),
        shared_proteins=sum(len(v) > 1 for v in owners.values()), reference_blocks=len(groups),
        merged_label_groups=sorted(group for group in groups if len(group) > 1))


def aggregate(raw, mapping):
    cells, statuses = defaultdict(Counter), defaultdict(set)
    with gzip.open(raw, "rt", newline="") as stream:
        for row in csv.reader(stream, delimiter="\t"):
            require(len(row) == 7, "Expected seven raw columns")
            a, b, category, fa, fb, _, _ = row
            require(a != b and category in CATEGORIES and fa in mapping and fb in mapping, "Invalid raw row")
            pair = tuple(sorted((a, b)))
            require(category not in statuses[pair], "Duplicate category/pair")
            statuses[pair].add(category)
            block_pair = tuple(sorted((mapping[fa], mapping[fb])))
            require(block_pair[0] == block_pair[1] or category == "FP", "Truth crosses reference blocks")
            cells[block_pair][category] += 1
    return cells, statuses


def ratios(counts):
    require(set(counts) == set(CATEGORIES) and all(type(v) is int and v >= 0 for v in counts.values()),
        "Invalid category counts")
    tp, fp, fn = (counts[c] for c in CATEGORIES)
    require(tp + fp > 0 and tp + fn > 0, "Undefined pooled ratio")
    return dict(precision=tp / (tp + fp), recall=tp / (tp + fn), f1=2 * tp / (2 * tp + fp + fn))


def close(actual, expected):
    require(set(actual) == set(expected) and all(type(expected[k]) in (int, float)
        and math.isfinite(expected[k]) and math.isclose(actual[k], expected[k], rel_tol=0, abs_tol=1e-12)
        for k in actual), "Pooled native metric mismatch")


def write_table(path, header, rows):
    with path.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def status(categories):
    return "+".join(c for c in CATEGORIES if c in categories) or "not_scored"


def export(snapshot_path, snapshot_sha, historical_path, historical_sha, output):
    output = Path(output).absolute()
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    checked = [record(__file__), record(PROTOCOL), record(Path(__file__).with_name("audit_qfo_vgnc_mapping.py")),
        record(Path(__file__).with_name("prepare_ob_candidate_neighborhood.py"))]
    snapshot_ref, historical_ref = record(snapshot_path), record(historical_path)
    require(snapshot_ref["sha256"] == snapshot_sha and historical_ref["sha256"] == historical_sha,
        "Changed frozen report")
    snapshot, historical = load(snapshot_ref, checked), load(historical_ref, checked)
    require(snapshot.get("schema") == "native_qfo_scientific_reporting_snapshot_v1"
        and snapshot.get("new_scoring_or_admission") is False and snapshot.get("publication_ready") is False
        and snapshot.get("recovered_inference_resources_admitted") is False, "Unexpected scientific snapshot scope")
    check(snapshot["source"])
    checked.append(snapshot["source"])
    load(snapshot["plan"], checked)
    require(historical.get("status") == "corrected_vgnc_native_rows_mapped_to_reference_blocks"
        and historical.get("uncertainty_admitted") is False, "Unexpected historical block scope")
    reference = unique_record(historical["checked_records"], "/vgnc-orthologs.txt.gz")
    check(reference)
    checked.append(reference)
    check(historical["reference_table"])
    checked.append(historical["reference_table"])
    truth, labels = reference_data(Path(reference["path"]))
    mapping, summary = reference_blocks(truth)
    require(summary == historical["reference"] and len(truth) == historical["reference_pairs"],
        "Changed reference-block reconstruction")
    admitted = [r for r in snapshot["rows"] if r.get("accuracy_admitted") is True]
    require([(r["index"], r["cell"]) for r in admitted] == list(zip((6, 7), CELLS)),
        "Require exactly the two completed native cells")
    output.mkdir(parents=True)
    methods, all_statuses, common_annotations = [], [], None
    for row in admitted:
        require(row["admission"] in snapshot["evidence"], "Snapshot admission not bound")
        admission = load(row["admission"], checked)
        require(admission.get("accuracy_admitted") is True and admission.get("publication_ready") is False
            and (admission["native_index"], admission["cell"], admission["participant"])
                == (row["index"], row["cell"], row["participant"]), "Native accuracy identity differs")
        recovered = row["cell"] == CELLS[1]
        require(admission["schema"] == ("measurement_failed_native_qfo_admission_v1" if recovered
            else "full_native_factorial_qfo_admission_v1") and admission["status"] == (
            "measurement_failed_native_qfo_assessment_admitted" if recovered
            else "full_native_factorial_qfo_assessment_admitted"), "Native admission schema/status differs")
        if recovered:
            require(row.get("timing_eligible") is False and row.get("timing_admitted") is False
                and admission.get("resources", "missing") is None
                and all(admission.get(k) is False for k in ("scientific_timings_admitted",
                    "eligible_for_timing_comparison", "original_native_scheduler_success")), "Failed timing relabeled")
        require(admission["execution_report"] in admission["checked_records"], "Execution not inventoried")
        execution = load(admission["execution_report"], checked)
        require(execution.get("exit_code") == 0 and execution.get("status") == "process_succeeded_pending_independent_admission"
            and (execution["native_index"], execution["cell"]) == (row["index"], row["cell"]),
            "Native assessment execution differs")
        participant = row["participant"]
        refs = [unique_record(admission["checked_records"], suffix) for suffix in (
            "/other/" + participant + ".db", "/results/VGNC/VGNC.json",
            "/VGNC_" + participant.replace("_", "-").replace(" ", "-") + "_raw.txt.gz")]
        require(all(r in execution["outputs"] for r in refs), "Artifact not originally inventoried by execution")
        for ref in refs:
            check(ref)
            checked.append(ref)
        database, metric, raw = (Path(r["path"]) for r in refs)
        pairs, annotations, digest, aliases = mapped_reference(database, truth, labels)
        require(digest == historical["reference_mapping_sha256"]
            and (common_annotations is None or annotations == common_annotations), "Reference mapping differs")
        common_annotations = annotations
        validation = validate_raw(raw, pairs, annotations)
        cells, statuses = aggregate(raw, mapping)
        totals = {c: sum(v[c] for v in cells.values()) for c in CATEGORIES}
        require(totals == validation["counts"], "Block count mismatch")
        observed = ratios(totals)
        inline = json.loads(metric.read_text())["datalink"]["inline_data"]
        native = inline["challenge_participants"]
        require(len(native) == 1 and native[0]["participant_id"] == participant
            and (inline["visualization"]["x_axis"], inline["visualization"]["y_axis"]) == ("TPR", "PPV"),
            "Native participant or axes differ")
        native = native[0]
        detail = row["endpoint_details"]["VGNC"]
        close(observed, dict(precision=detail["precision"], recall=detail["recall"], f1=row["scores"]["VGNC"]))
        close({k: observed[k] for k in ("precision", "recall")},
            dict(precision=native["metric_y"], recall=native["metric_x"]))
        require(native == admission["assessment"]["endpoints"]["VGNC"]["native_participant"],
            "Aggregate differs from original assessment")
        table = output / (row["cell"] + ".tsv")
        write_table(table, ["block_left", "block_right", *CATEGORIES],
            [[a, b, *[v[c] for c in CATEGORIES]] for (a, b), v in sorted(cells.items())])
        methods.append(dict(cell=row["cell"], index=row["index"], participant=participant,
            prediction_semantics=row["prediction_semantics"], measurement_status=row["measurement_status"],
            admission=row["admission"], execution=admission["execution_report"],
            database=refs[0], aggregate=refs[1], raw=refs[2], table=record(table), counts=totals,
            metrics=observed, validation=validation, alias_rows=aliases,
            selected_reference_mapping_sha256=digest, full_database_hash_checked=True,
            prediction_edges_requeried=False, nonzero_cells=len(cells),
            nonzero_cross_block_cells=sum(a != b for a, b in cells),
            within_block_false_positives=sum(v["FP"] for (a, b), v in cells.items() if a == b),
            cross_block_false_positives=sum(v["FP"] for (a, b), v in cells.items() if a != b)))
        all_statuses.append(statuses)
    proteins, taxa, asserted = defaultdict(set), defaultdict(set), Counter()
    for accession, (label, taxon) in common_annotations.items():
        proteins[mapping[label]].add(accession)
        taxa[mapping[label]].add(taxon)
    for family in truth.values():
        asserted[mapping[family]] += 1
    inventory = output / "reference_blocks.tsv"
    write_table(inventory, ["block", "proteins", "asserted_pairs", "species"],
        [[block, len(proteins[block]), asserted[block], ",".join(sorted(taxa[block]))]
            for block in sorted(set(mapping.values()))])
    require(inventory.read_bytes() == Path(historical["reference_table"]["path"]).read_bytes(),
        "Full reference inventory differs")
    transition_counts, transitions = Counter(), []
    for a, b in sorted(set(all_statuses[0]) | set(all_statuses[1])):
        left, right = (status(s.get((a, b), set())) for s in all_statuses)
        blocks = sorted(mapping[common_annotations[g][0]] for g in (a, b))
        transitions.append([a, b, *blocks, left, right])
        transition_counts[left, right] += 1
    transition_table = output / "pair_transitions.tsv"
    write_table(transition_table, ["protein_left", "protein_right", "block_left", "block_right", *CELLS], transitions)
    for ref in checked:
        check(ref)
    result = dict(schema="native_qfo_vgnc_blocks_v1", status="native_scored_rows_decomposed",
        source=record(__file__), protocol=record(PROTOCOL), snapshot=snapshot_ref, historical_blocks=historical_ref,
        checked_records=checked, reference_record=reference, reference=summary, reference_pairs=len(truth),
        reference_table=record(inventory), methods=methods, transition_table=record(transition_table),
        union_scored_pairs=len(transitions), transition_counts=[dict(r0=a, r1=b, pairs=n)
            for (a, b), n in sorted(transition_counts.items())],
        differences={k: methods[1]["metrics"][k] - methods[0]["metrics"][k] for k in methods[0]["metrics"]},
        uncertainty_admitted=False, new_scoring_or_admission=False, publication_ready=False,
        failed_r1_timing_remains_ineligible=True,
        zero_cell_convention="Omitted fixed unordered block cells contribute zero counts, not proven biological ineligibility.",
        limitations=["Development-exposed native cells, not historical selected-default configurations.",
            "Current-byte checks against original inventories do not prove uninterrupted artifact integrity.",
            "Native last-rowid alias reconstruction is tested here; native SQL has no general alias-order guarantee.",
            "Prediction edges were not queried; completeness of omitted false positives is not independently established.",
            "Category overlap is retained; not_scored is not a true-negative or biological non-orthology label.",
            "Reference overlap blocks do not establish independent biological sampling units or valid confidence intervals.",
            "No causal evolutionary claim, new method/default, total-HMM ablation or recovered timing admission.",
            "Diagnostic timing is shared-host postprocessing, not native inference timing; contention impact is unknown and potentially tool-dependent."])
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent.parent
    parser.add_argument("--snapshot", type=Path, default=root / "benchmark_tools/results/native_qfo_scientific_scores_20261006_v1/report.json")
    parser.add_argument("--snapshot-sha256", default=SNAPSHOT_SHA)
    parser.add_argument("--historical", type=Path, default=root / "benchmark_tools/results/corrected_vgnc_blocks_20260926.json")
    parser.add_argument("--historical-sha256", default=HISTORICAL_SHA)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = export(args.snapshot, args.snapshot_sha256, args.historical, args.historical_sha256, args.output)
    print(json.dumps(dict(status=result["status"], methods=len(result["methods"]), union_scored_pairs=result["union_scored_pairs"])))

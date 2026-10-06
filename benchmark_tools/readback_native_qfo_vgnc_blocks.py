"""Independent stdlib VGNC raw/SQLite/union-find/ratio readback."""

import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3

CATEGORIES = ("TP", "FP", "FN")
CELLS = ("p0_c0_r0", "p0_c0_r1")


def need(condition, message):
    if not condition:
        raise ValueError(message)


def identity(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def verify(ref):
    need(identity(ref["path"]) == ref, "Changed evidence: " + ref["path"])


def read_json(ref):
    verify(ref)
    return json.loads(Path(ref["path"]).read_text())


def reconstruct(reference):
    assertions, labels, owners, parent = {}, {}, defaultdict(set), {}
    with gzip.open(reference, "rt", newline="") as stream:
        for row in csv.reader(stream, delimiter="\t"):
            need(len(row) == 3, "Bad reference columns")
            a, b, family = int(row[0]), int(row[1]), row[2]
            pair = (min(a, b), max(a, b))
            need(a != b and pair not in assertions, "Duplicate/self reference pair")
            assertions[pair] = family
            parent.setdefault(family, family)
            for protein in pair:
                owners[protein].add(family)
                labels[protein] = family

    def find(label):
        while label != parent[label]:
            parent[label] = parent[parent[label]]
            label = parent[label]
        return label

    for families in owners.values():
        ordered = sorted(families)
        for family in ordered[1:]:
            a, b = find(ordered[0]), find(family)
            parent[max(a, b)] = min(a, b)
    mapping = {family: find(family) for family in parent}
    groups = defaultdict(list)
    for family, block in sorted(mapping.items()):
        groups[block].append(family)
    summary = dict(family_labels=len(parent), reference_proteins=len(owners),
        shared_proteins=sum(len(v) > 1 for v in owners.values()), reference_blocks=len(groups),
        merged_label_groups=sorted(v for v in groups.values() if len(v) > 1))
    return assertions, labels, mapping, summary


def mapped_database(path, labels):
    selected, accession_owner = {}, {}
    extras = dict(identical=0, alias=0)
    with sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True) as connection:
        cursor = connection.execute("SELECT rowid, prot_nr, uniprot_id, species FROM proteomes ORDER BY rowid")
        for _, protein, accession, species in cursor:
            if protein not in labels:
                continue
            value = [accession, labels[protein], species]
            if selected.get(protein) == value:
                extras["identical"] += 1
                continue
            need(accession not in accession_owner or accession_owner[accession] == protein,
                "Non-bijective accession")
            if protein in selected:
                need(selected[protein][2] == species, "Conflicting species")
                extras["alias"] += 1
            accession_owner[accession] = protein
            selected[protein] = value
    need(set(selected) == set(labels), "Incomplete reference mapping")
    digest = hashlib.sha256(json.dumps(sorted(selected.items()), separators=(",", ":")).encode()).hexdigest()
    return selected, digest, extras


def raw_counts(path, asserted, annotations, mapping):
    categories = {c: set() for c in CATEGORIES}
    cells, pairs = defaultdict(Counter), defaultdict(set)
    species_families = defaultdict(set)
    for label, species in annotations.values():
        species_families[species].add(label)
    rows = 0
    with gzip.open(path, "rt") as stream:
        for line in stream:
            row = line.rstrip("\n").split("\t")
            need(len(row) == 7, "Bad raw columns")
            a, b, category, fa, fb, sa, sb = row
            need(category in categories and a != b, "Invalid category/self pair")
            need(annotations.get(a) == (fa, sa) and annotations.get(b) == (fb, sb), "Raw annotation mismatch")
            pair = tuple(sorted((a, b)))
            need(pair not in categories[category], "Duplicate category/pair")
            if category != "FP":
                need(pair in asserted, "Unasserted TP/FN")
            else:
                need(fa != fb and fa in species_families[sb] and fb in species_families[sa], "Ineligible FP")
            categories[category].add(pair)
            pairs[pair].add(category)
            blocks = tuple(sorted((mapping[fa], mapping[fb])))
            need(blocks[0] == blocks[1] or category == "FP", "Cross-block truth")
            cells[blocks][category] += 1
            rows += 1
    need(not categories["TP"] & categories["FN"] and categories["TP"] | categories["FN"] == asserted,
        "Incomplete truth partition")
    validation = dict(counts={c: len(categories[c]) for c in CATEGORIES},
        tp_fp_overlap=len(categories["TP"] & categories["FP"]), fn_fp_overlap=len(categories["FN"] & categories["FP"]))
    return cells, pairs, validation, rows


def table(ref, header, expected):
    verify(ref)
    with Path(ref["path"]).open(newline="") as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    need(rows and rows[0] == header and rows[1:] == [[str(v) for v in row] for row in expected],
        "Complete table mismatch: " + ref["path"])


def near(actual, expected):
    need(type(expected) in (int, float) and abs(float(actual) - expected) <= 1e-12, "Ratio mismatch")


def review(path, digest):
    report_ref = identity(path)
    need(report_ref["sha256"] == digest, "Changed report")
    report = read_json(report_ref)
    need(report.get("schema") == "native_qfo_vgnc_blocks_v1" and report.get("status") == "native_scored_rows_decomposed"
        and all(report.get(k) is False for k in ("uncertainty_admitted", "new_scoring_or_admission", "publication_ready"))
        and report.get("failed_r1_timing_remains_ineligible") is True, "Changed diagnostic scope")
    checked = report["checked_records"]
    need(report["source"] == identity(Path(__file__).with_name("export_native_qfo_vgnc_blocks.py")),
        "Exporter source changed")
    for ref in checked:
        verify(ref)
    need(all(ref in checked for ref in (report["source"], report["protocol"], report["snapshot"],
        report["historical_blocks"], report["reference_record"])), "Unbound inputs")
    snapshot, historical = read_json(report["snapshot"]), read_json(report["historical_blocks"])
    need(snapshot["schema"] == "native_qfo_scientific_reporting_snapshot_v1"
        and snapshot["new_scoring_or_admission"] is False and snapshot["publication_ready"] is False
        and snapshot["recovered_inference_resources_admitted"] is False, "Changed snapshot scope")
    need(snapshot["source"] in checked and snapshot["plan"] in checked
        and historical["reference_table"] in checked and report["reference_record"] in historical["checked_records"],
        "Missing direct source/reference binding")
    assertions, labels, mapping, summary = reconstruct(report["reference_record"]["path"])
    need(summary == historical["reference"] == report["reference"]
        and len(assertions) == historical["reference_pairs"] == report["reference_pairs"], "Reference mismatch")
    methods = report["methods"]
    need([m["cell"] for m in methods] == list(CELLS), "Wrong native cells")
    snapshot_rows = [r for r in snapshot["rows"] if r.get("accuracy_admitted") is True]
    need([(r["index"], r["cell"]) for r in snapshot_rows] == list(zip((6, 7), CELLS)), "Changed admitted cohort")
    statuses, reference_annotations, metrics, raw_rows = [], None, [], 0
    for method, row in zip(methods, snapshot_rows):
        need(all(method[k] == row[k] for k in ("index", "cell", "participant", "prediction_semantics", "measurement_status", "admission"))
            and row["admission"] in snapshot["evidence"] and method["admission"] in checked, "Admission row mismatch")
        admission = read_json(method["admission"])
        need(admission["accuracy_admitted"] is True and admission["publication_ready"] is False
            and admission["native_index"] == method["index"] and admission["cell"] == method["cell"]
            and admission["participant"] == method["participant"] and admission["execution_report"] == method["execution"]
            and method["execution"] in checked and method["execution"] in admission["checked_records"], "Changed admission identity")
        if method["cell"] == CELLS[1]:
            need(row.get("timing_eligible") is False and row.get("timing_admitted") is False
                and admission.get("resources", "missing") is None
                and all(admission.get(k) is False for k in ("scientific_timings_admitted", "eligible_for_timing_comparison",
                    "original_native_scheduler_success")), "Failed timing changed")
        execution = read_json(method["execution"])
        need(execution["exit_code"] == 0 and execution["native_index"] == method["index"]
            and execution["cell"] == method["cell"], "Execution differs")
        for key in ("database", "raw", "aggregate"):
            need(method[key] in checked and method[key] in admission["checked_records"]
                and method[key] in execution["outputs"], "Original artifact inventory missing")
        selected, mapping_digest, aliases = mapped_database(method["database"]["path"], labels)
        need(mapping_digest == historical["reference_mapping_sha256"] == method["selected_reference_mapping_sha256"]
            and aliases == method["alias_rows"] and method["full_database_hash_checked"] is True
            and method["prediction_edges_requeried"] is False, "Mapping identity/scope mismatch")
        annotations = {v[0]: (v[1], v[2]) for v in selected.values()}
        need(reference_annotations is None or annotations == reference_annotations, "Cross-cell mapping differs")
        reference_annotations = annotations
        asserted = {tuple(sorted((selected[a][0], selected[b][0]))) for a, b in assertions}
        cells, pairs, validation, count = raw_counts(method["raw"]["path"], asserted, annotations, mapping)
        raw_rows += count
        need(validation == method["validation"] and validation["counts"] == method["counts"], "Counts differ")
        table(method["table"], ["block_left", "block_right", *CATEGORIES],
            [[a, b, *[v[c] for c in CATEGORIES]] for (a, b), v in sorted(cells.items())])
        need(method["nonzero_cells"] == len(cells)
            and method["nonzero_cross_block_cells"] == sum(a != b for a, b in cells)
            and method["within_block_false_positives"] == sum(v["FP"] for (a, b), v in cells.items() if a == b)
            and method["cross_block_false_positives"] == sum(v["FP"] for (a, b), v in cells.items() if a != b),
            "Block summaries differ")
        tp, fp, fn = (validation["counts"][c] for c in CATEGORIES)
        observed = dict(precision=Fraction(tp, tp + fp), recall=Fraction(tp, tp + fn), f1=Fraction(2 * tp, 2 * tp + fp + fn))
        for key, value in observed.items():
            near(value, method["metrics"][key])
            near(value, row["scores"]["VGNC"] if key == "f1" else row["endpoint_details"]["VGNC"][key])
        inline = read_json(method["aggregate"])["datalink"]["inline_data"]
        native = inline["challenge_participants"]
        need(len(native) == 1 and native[0]["participant_id"] == method["participant"]
            and inline["visualization"]["x_axis"] == "TPR" and inline["visualization"]["y_axis"] == "PPV"
            and native[0] == admission["assessment"]["endpoints"]["VGNC"]["native_participant"], "Native aggregate differs")
        near(observed["precision"], native[0]["metric_y"])
        near(observed["recall"], native[0]["metric_x"])
        metrics.append(observed)
        statuses.append(pairs)
    proteins, taxa, asserted_counts = defaultdict(set), defaultdict(set), Counter()
    for accession, (family, species) in reference_annotations.items():
        proteins[mapping[family]].add(accession)
        taxa[mapping[family]].add(species)
    for family in assertions.values():
        asserted_counts[mapping[family]] += 1
    reference_rows = [[b, len(proteins[b]), asserted_counts[b], ",".join(sorted(taxa[b]))] for b in sorted(set(mapping.values()))]
    for ref in (report["reference_table"], historical["reference_table"]):
        table(ref, ["block", "proteins", "asserted_pairs", "species"], reference_rows)
    transitions, counts = [], Counter()
    for a, b in sorted(set(statuses[0]) | set(statuses[1])):
        values = ["+".join(c for c in CATEGORIES if c in p.get((a, b), set())) or "not_scored" for p in statuses]
        blocks = sorted(mapping[reference_annotations[g][0]] for g in (a, b))
        transitions.append([a, b, *blocks, *values])
        counts[tuple(values)] += 1
    table(report["transition_table"], ["protein_left", "protein_right", "block_left", "block_right", *CELLS], transitions)
    need(len(transitions) == report["union_scored_pairs"] and report["transition_counts"]
        == [dict(r0=a, r1=b, pairs=n) for (a, b), n in sorted(counts.items())], "Transition counts differ")
    for key in metrics[0]:
        near(metrics[1][key] - metrics[0][key], report["differences"][key])
    for ref in checked:
        verify(ref)
    verify(report_ref)
    return dict(schema="native_qfo_vgnc_blocks_readback_v1", status="complete_native_decomposition_independently_verified",
        source=identity(__file__), report=report_ref, checked_input_records=len(checked), raw_rows=raw_rows,
        reference_rows=len(reference_rows), transition_pairs=len(transitions), transitions=report["transition_counts"],
        methods=[dict(cell=m["cell"], counts=m["counts"], sparse_cells=m["nonzero_cells"]) for m in methods],
        exporter_or_reference_helpers_imported=False, uncertainty_admitted=False, new_scoring_or_admission=False,
        failed_r1_timing_remains_ineligible=True, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = review(args.report, args.report_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write("\n")
    print(json.dumps(dict(status=result["status"], raw_rows=result["raw_rows"], transition_pairs=result["transition_pairs"])))

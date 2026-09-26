"""Map all retained corrected VGNC scored rows to fixed reference blocks."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path

from benchmark_tools.audit_qfo_vgnc_mapping import reference_data, mapped_reference, validate_raw
from benchmark_tools.audit_vgnc_family_dependencies import blocks, REFERENCE_SHA
from benchmark_tools.diagnose_vgnc_block_influence import metrics
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

MANIFEST = "benchmark_tools/results/qfo_corrected_comparison_20260926_v7/manifest.json"
MANIFEST_SHA = "042aa221d01554114a1b0b413ca3e2ff56dad5f8c06362c69a1bf49f2de642fc"
CATEGORIES = ("TP", "FP", "FN")


def aggregate(rows, mapping):
    cells, seen = defaultdict(Counter), set()
    for row in rows:
        if len(row) != 7:
            raise ValueError("Expected seven raw columns")
        a, b, category, fa, fb, sa, sb = row
        if a == b or category not in CATEGORIES or fa not in mapping or fb not in mapping:
            raise ValueError("Invalid scored row")
        key = (category, *sorted((a, b)))
        if key in seen:
            raise ValueError("Duplicate category/pair")
        seen.add(key)
        left, right = sorted((mapping[fa], mapping[fb]))
        if left != right and category != "FP":
            raise ValueError("Truth crosses reference blocks")
        cells[left, right][category] += 1
    return cells


def write_cells(path, cells):
    with path.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["block_left", "block_right", *CATEGORIES])
        for (left, right), counts in sorted(cells.items()):
            writer.writerow([left, right, *[counts[c] for c in CATEGORIES]])


def run(repo, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    manifest_record = record(repo / MANIFEST)
    if manifest_record["sha256"] != MANIFEST_SHA:
        raise ValueError("Changed frozen comparator manifest")
    manifest = json.loads((repo / MANIFEST).read_text())
    if len(manifest["methods"]) != 8 or any(m["status"] != "admitted" for m in manifest["methods"]):
        raise ValueError("Expected eight admitted methods")
    reference = repo / "qfo_benchmark/benchmark-webservice/reference_data/2020/vgnc-orthologs.txt.gz"
    reference_record = record(reference)
    if reference_record["sha256"] != REFERENCE_SHA:
        raise ValueError("Changed VGNC reference")
    truth, labels = reference_data(reference)
    mapping, summary = blocks(truth)
    checked = [manifest_record, reference_record, record(__file__)]
    checked.extend(record(Path(__file__).with_name(name)) for name in (
        "audit_qfo_vgnc_mapping.py", "audit_vgnc_family_dependencies.py", "diagnose_vgnc_block_influence.py"))
    output.mkdir(parents=True)
    methods, mapping_digest, common_annotations = [], None, None
    for method in manifest["methods"]:
        check(method["admission"])
        checked.append(method["admission"])
        admission = json.loads(Path(method["admission"]["path"]).read_text())
        candidates = [r for r in admission["metric_files"] if str(r["path"]).endswith("/results/VGNC/VGNC.json")]
        if len(candidates) != 1:
            raise ValueError("Ambiguous native VGNC aggregation")
        metric_record = candidates[0]
        check(metric_record)
        checked.append(metric_record)
        metric_path = Path(metric_record["path"])
        participant_rows = json.loads(metric_path.read_text())["datalink"]["inline_data"]["challenge_participants"]
        if len(participant_rows) != 1:
            raise ValueError("Ambiguous native participant")
        participant = participant_rows[0]["participant_id"]
        raw = metric_path.parent / ("VGNC_" + participant.replace(" ", "-").replace("_", "-") + "_raw.txt.gz")
        database = metric_path.parents[2] / "other" / (participant + ".db")
        raw_record = record(raw)
        checked.append(raw_record)
        before = database.stat()
        pairs, annotations, digest, aliases = mapped_reference(database, truth, labels)
        after = database.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError("Prediction database changed during reference readback")
        if mapping_digest is not None and (digest != mapping_digest or annotations != common_annotations):
            raise ValueError("Reference mappings differ across methods")
        mapping_digest, common_annotations = digest, annotations
        validation = validate_raw(raw, pairs, annotations)
        with gzip.open(raw, "rt", newline="") as stream:
            cells = aggregate(csv.reader(stream, delimiter="\t"), mapping)
        totals = {c: sum(row[c] for row in cells.values()) for c in CATEGORIES}
        if totals != validation["counts"]:
            raise ValueError("Block aggregation lost native scored rows")
        observed = metrics(totals)
        expected = dict(precision=method["details"]["VGNC"]["precision"],
                        recall=method["details"]["VGNC"]["recall"], f1=method["scores"]["VGNC"])
        if any(abs(observed[c] - expected[c]) > 1e-12 for c in expected):
            raise ValueError("Raw/block metrics differ from admitted scores")
        table = output / (method["key"] + ".tsv")
        write_cells(table, cells)
        cross = [(pair, count) for pair, count in cells.items() if pair[0] != pair[1]]
        methods.append(dict(key=method["key"], label=method["label"], participant=participant,
            prediction_semantics=method["prediction_semantics"], raw=raw_record, table=record(table),
            counts=totals, metrics=observed, raw_validation=validation,
            cross_block_false_positives=sum(v["FP"] for _, v in cross),
            within_block_false_positives=sum(v["FP"] for pair, v in cells.items() if pair[0] == pair[1]),
            nonzero_cross_block_cells=len(cross), nonzero_cells=len(cells),
            database=dict(path=str(database), bytes=before.st_size, mtime_ns=before.st_mtime_ns,
                          selected_reference_mapping_sha256=digest, alias_rows=aliases,
                          full_database_hash_checked=False, prediction_edges_requeried=False)))
    proteins, species, reference_counts = defaultdict(set), defaultdict(set), Counter()
    for accession, (label, taxon) in common_annotations.items():
        proteins[mapping[label]].add(accession)
        species[mapping[label]].add(taxon)
    for _, family in truth.items():
        reference_counts[mapping[family]] += 1
    inventory = output / "reference_blocks.tsv"
    with inventory.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["block", "proteins", "asserted_pairs", "species"])
        for block in sorted(set(mapping.values())):
            writer.writerow([block, len(proteins[block]), reference_counts[block], ",".join(sorted(species[block]))])
    for item in checked:
        check(item)
    result = dict(status="corrected_vgnc_native_rows_mapped_to_reference_blocks", checked_records=checked,
        reference=summary, reference_pairs=len(truth), reference_table=record(inventory), methods=methods,
        reference_mapping_sha256=mapping_digest,
        protein_count_histogram=dict(sorted(Counter(map(len, proteins.values())).items())),
        reference_pair_count_histogram=dict(sorted(Counter(reference_counts.values()).items())),
        zero_cell_convention="All fixed reference-block diagonal and unordered off-diagonal cells exist; omitted cells have zero TP/FP/FN contributions, not necessarily equal biological eligibility.",
        uncertainty_admitted=False, publication_ready=False,
        limitations=["Raw score decomposition, not a new prediction-database rescore or historical execution proof.",
            "Database checks cover selected reference mappings and size/mtime stability, not full content hashes.",
            "Raw hashes are established in this readback; score admission bound aggregates, not these supplemental raw tables.",
            "Reference-defined overlap blocks do not establish independent evolutionary or prediction units.",
            "No synthetic variance estimator is applied and no biological confidence interval is admitted."])
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.output.absolute())

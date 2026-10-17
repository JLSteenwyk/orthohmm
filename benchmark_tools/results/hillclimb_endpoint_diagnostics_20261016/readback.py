"""Verify endpoint diagnostics directly against frozen JSON locations and TSVs."""

import collections
import hashlib
import json
from pathlib import Path
import csv
import sys
import xml.etree.ElementTree as ET


def pin(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            digest.update(block)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def follow(value, path):
    for key in path:
        value = value[key]
    return value


def main(root):
    report_path = root / "report.json"
    report = json.loads(report_path.read_text())
    pins = [*report["inputs"], *report["tables"]]
    for item in pins:
        assert pin(item["path"]) == item
    sources = {item["path"]: json.loads(Path(item["path"]).read_text()) for item in report["inputs"][:2]}
    summary = sources[report["inputs"][0]["path"]]
    methods = {m["key"]: m for m in summary["methods"]}
    admitted = {"retained_registered_result", "independent_content_readback_verified"}
    observed, latest = collections.Counter(), {}
    largest_harmonic_error = 0
    for row in report["history"]:
        method = methods[row["method"]]
        cell = method["cells"][row["cell_index"]]
        assert row["cell_status"] == cell["status"] and row["dataset"] == cell["dataset"]
        observed[(row["method"], row["cell_index"])] += 1
        for key in ("genes", "species", "groups", "pairs"):
            assert row["native_" + key] == cell.get("coverage", {}).get(key)
        if row["cell_status"] not in admitted:
            assert row["score"] is None
            assert row.get("precision") is None and row.get("recall") is None
            continue
        score = follow(sources[row["score_source_document"]], row["score_source_location"])
        assert score == row["score"] == method["scores"][row["dataset"] + ": " + row["endpoint"]]
        entry = follow(sources[row["detail_source_document"]], row["detail_source_location"])
        if row["dataset"] == "OrthoBench" and row["cell_status"] == "retained_registered_result":
            for field in ("precision", "recall"):
                assert row[field] == entry[field + "_percent"] * 0.01 or abs(row[field] - entry[field + "_percent"] * 0.01) < 1e-15
            assert row["benchmark_genes_in_orthogroups"] == entry["genes_in_orthogroups"]
            assert abs(row["benchmark_reported_gene_coverage"] - entry["coverage_percent"] / 100) < 1e-15
            assert row["exact_refogs"] == entry["exact_refogs"]
            assert abs(entry["f_score_percent"] / 100 - score) < 1e-12
            retained_row = follow(sources[row["detail_source_document"]], row["detail_source_location"][:-1])
            prediction = entry.get("prediction_provenance") or retained_row["orthobench_supplemental_readback"]["prediction"]
            assert prediction == cell["retained_details"]["retained_ob_metadata"]["prediction"]
        else:
            assert row.get("precision") == entry.get("precision")
            assert row.get("recall") == entry.get("recall")
        if row["dataset"] == "QfO":
            assert row.get("assessed_relations") == entry.get("assessed_relations")
            retained_details = cell.get("retained_details", {})
            for field in ("submitted_pairs", "retained_pairs"):
                assert row.get(field) == retained_details.get(field)
        elif row["dataset"] == "ThreeKingdoms":
            for field in ("reference_gene_coverage", "reference_genes", "reference_genes_in_prediction",
                          "true_positive_gene_pairs", "false_positive_gene_pairs", "false_negative_gene_pairs"):
                assert row[field] == entry[field]
            assert row["benchmark_predicted_orthogroups"] == entry["predicted_orthogroups"]
            assert row["score"] == entry["f_score"]
        elif row["cell_status"] != "retained_registered_result":
            assert row.get("exact_refogs") == entry.get("family_details", {}).get("exact_refogs")
        if row["statistic"] == "conditional_similarity":
            assert row.get("precision") is None and row.get("recall") is None
            if "similarity" in entry:
                assert entry["similarity"] == row["score"]
        else:
            p, r = row["precision"], row["recall"]
            harmonic = 2 * p * r / (p + r) if p + r else 0
            largest_harmonic_error = max(largest_harmonic_error, abs(harmonic - row["score"]))
            assert abs(harmonic - row["score"]) < 1e-10
        latest[(row["method"], row["dataset"], row["endpoint"])] = row
    all_cells = {(m["key"], i): 6 if c["dataset"] == "QfO" else 1
                 for m in summary["methods"] for i, c in enumerate(m["cells"])}
    assert dict(observed) == all_cells
    assert len(all_cells) == 91 and len(methods) == 28
    current_index = {}
    for row in report["current"]:
        key = row["method"], row["dataset"], row["endpoint"]
        assert key not in current_index
        current_index[key] = row
        if key in latest:
            assert row == latest[key]
        else:
            assert row["cell_status"] == "no_admitted_endpoint" and row["score"] is None
            assert row["dataset"] + ": " + row["endpoint"] not in methods[row["method"]]["scores"]
    assert len(current_index) == 224
    assert sum(row.get("score") is not None for row in report["current"]) == 213
    for item, rows in zip(report["tables"], (report["current"], report["history"])):
        with Path(item["path"]).open() as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            table = list(reader)
            columns = reader.fieldnames
        assert len(table) == len(rows)
        for table_row, json_row in zip(table, rows):
            for key in columns:
                value = json_row.get(key)
                text = "NA" if value is None else json.dumps(value) if isinstance(value, list) else str(value)
                assert table_row[key] == text
    tests = Path("benchmarks/results/hillclimb_diagnostic_export_tests_20261016.xml")
    suites = list(ET.parse(tests).getroot().iter("testsuite"))
    counts = {key: sum(int(s.attrib.get(key, 0)) for s in suites) for key in ("tests", "failures", "errors", "skipped")}
    assert counts == dict(tests=71, failures=0, errors=0, skipped=0)
    contrasts = []
    for row in report["current"]:
        if row["method"] != "canonical_collection:canonical":
            continue
        other = current_index[("orthofinder_3_1_5_full", row["dataset"], row["endpoint"])]
        contrasts.append(dict(dataset=row["dataset"], endpoint=row["endpoint"], statistic=row["statistic"],
            orthohmm=row["score"], orthofinder_full=other["score"],
            differences={key: row[key] - other[key] if row.get(key) is not None and other.get(key) is not None else None
                         for key in ("score", "precision", "recall")}))
    for item in pins:
        assert pin(item["path"]) == item
    receipt = dict(status="independent_endpoint_source_units_ownership_and_tsv_readback_verified",
        methods=len(methods), historical_cells=len(all_cells), current_rows=224, scored_current_rows=213,
        historical_endpoint_observations=len(report["history"]), missing_current_endpoints=11,
        scored_current_F1_rows=132, max_F1_harmonic_absolute_difference=largest_harmonic_error,
        canonical_vs_orthofinder_full=contrasts, checked=pins, report=pin(report_path), recipe=pin(__file__),
        tests={"counts": counts, "xml": pin(tests)}, defaults_promoted=False, goal_completed=False,
        limitations=["Independent direct-location readback and harmonic arithmetic check, not new raw biological scoring or independent validation.",
            "Historical precision/recall rounding is preserved; no F1 precision/recall inferred from F1 alone.",
            "All 91 historical cells retain their endpoint observations; current table uses the latest admitted cell and eleven missing endpoints remain NA.",
            "Canonical versus OrthoFinder point contrasts are descriptive. GO/EC/FAS populations differ; FAS is unseeded. No universal winner or statistical significance follows."])
    with (root / "independent_readback.json").open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"current_rows": 224, "historical_observations": len(report["history"]),
        "scored_current_rows": 213, "F1_rows": 132, "tests": counts}))


if __name__ == "__main__":
    main(Path(sys.argv[1]))

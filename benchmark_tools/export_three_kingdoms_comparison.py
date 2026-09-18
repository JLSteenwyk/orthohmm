"""Export audited BUSCO co-membership scores with explicit run provenance."""

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_three_kingdoms_pair_counts import count, compare
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.publication_comparison import METHODS
from benchmark_tools.run_simulation_methods import read_frozen

SOURCES = {
    "three_kingdoms_pair_counts_audit_20260918.json": "9a0ca8e1440acba1616d261fabd96c14801aecf3244ed6ed53bfc0135ab3bf1e",
    "three_kingdoms_sonic_matched_assessment_21796.json": "aeddd422a5c733ac40bef0c4953a6bda12ab2a751baeeae9e9c28cac472e9653",
    "three_kingdoms_method_inputs_20260918.json": "5dfd585193a6a576f2685388827bd6f8e3f5bad2ffcc715aaf665247cd2345dc",
}
SONIC = "sonicparanoid_2_0_9"


def assemble(historical, matched, inputs):
    if (historical["status"] != "historical_three_kingdoms_pair_counts_independently_verified"
            or matched["status"] != "matched_three_kingdoms_sonic_score_verified"
            or matched["accuracy_admitted"] is not True):
        raise ValueError("Require validated count reports")
    expected = {row[0] for row in METHODS}
    old = {row["key"]: row for row in historical["methods"]}
    provenance = {row["method"]: row for row in inputs["methods"]}
    if (set(old) != expected or set(provenance) != expected
            or len(historical["methods"]) != 8 or len(inputs["methods"]) != 8):
        raise ValueError("Missing or duplicate method")
    rows = []
    for key, label, _, _ in METHODS:
        contemporary = key == SONIC
        source = matched if contemporary else old[key]
        rows.append(dict(key=key, method=label, counts=source["counts"],
            groups=source["normalized"] if contemporary else source["groups"],
            run="contemporary matched input" if contemporary else "historical",
            input_status="matched input verified" if contemporary else
                provenance[key]["input_hash_manifest"] + "; historical consumption not proven",
            use="comparison", semantics="reference-gene group co-membership, including within-species pairs"))
    rows.append(dict(key=SONIC + "_historical", method="SonicParanoid 2.0.9 (historical)",
        counts=old[SONIC]["counts"], groups=old[SONIC]["groups"], run="historical",
        input_status="known Danio input mismatch", use="historical diagnostic only",
        semantics="reference-gene group co-membership, including within-species pairs"))
    for row in rows:
        c = row["counts"]
        if (c["reference_genes"] != 2035 or c["reference_orthogroups"] != 255
                or c["true_positive_gene_pairs"] + c["false_negative_gene_pairs"] != 7352):
            raise ValueError("Reference universe differs")
    return rows


def export(root, output):
    if output.exists():
        raise FileExistsError(output)
    base = root / "benchmark_tools/results"
    documents = [read_frozen(base / name, sha) for name, sha in SOURCES.items()]
    rows = assemble(*documents)
    evidence = [record(base / name) for name in SOURCES]
    historical, matched, _ = documents
    reference = next(item for item in historical["evidence"]
                     if Path(item["path"]).name == "reference_orthogroups.txt")
    if reference not in matched["checked_records"]:
        raise ValueError("Contemporary and historical reference identity differs")
    evidence.extend([reference, *historical["evidence"], *matched["checked_records"],
                     matched["score"], matched["normalized"], matched["conversion_log"], matched["source"]])
    for item in evidence:
        check(item)
    for row in rows:
        check(row["groups"])
        compare(count(reference["path"], row["groups"]["path"]), row["counts"])
    for item in evidence:
        check(item)
    result = dict(status="three_kingdoms_comparison_recounted", rows=rows, evidence=evidence,
        source=record(__file__), helpers=[record(Path(__file__).with_name(name)) for name in (
            "audit_three_kingdoms_pair_counts.py", "publication_comparison.py")],
        publication_ready=False, uniform_historical_input_consumption_proven=False,
        limitations=["Supplementary conserved-family endpoint, not genome-wide orthology accuracy.",
            "Predictions outside the reference universe are not penalized.",
            "Historical normalized-group count audits do not prove identical historical input consumption.",
            "Historical Sonic is diagnostic only; the new run does not isolate the cause of its difference.",
            "No confidence intervals or superiority tests are supplied by this descriptive table.",
            "OrthoFinder sequence-only is an MCL checkpoint; FastOMA uses a supplied OrthoFinder tree."])
    output.mkdir(parents=True, exist_ok=False)
    (output / "comparison.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    fields = ["method", "run", "use", "precision", "recall", "f_score", "reference_gene_coverage", "input_status"]
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            flat = {**row, **row["counts"]}
            writer.writerow({key: flat[key] for key in fields})
    lines = ["# Three Kingdoms: Audited Group Co-membership", "",
             "| Method | Run | Precision | Recall | F1 | Reference coverage |",
             "|---|---|---:|---:|---:|---:|"]
    for row in rows:
        c = row["counts"]
        values = " | ".join(f"{c[k]:.6f}" for k in ("precision", "recall", "f_score", "reference_gene_coverage"))
        lines.append(f"| {row['method']} | {row['run']} | {values} |")
    lines.extend(["", "Historical Sonic is diagnostic only because of its known input mismatch.", ""])
    lines.extend("- " + limit for limit in result["limitations"])
    lines.extend(["", "## Input Evidence", ""])
    lines.extend(f"- {row['method']}: {row['input_status']}." for row in rows)
    (output / "scores.md").write_text("\n".join(lines) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    export(args.root.resolve(), args.output.resolve())

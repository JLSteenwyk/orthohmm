"""Project admitted native family counts into independently verified distance bins."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0")
METRICS = ("F1", "PPV", "TPR")
CONTRASTS = (("R_at_P0_C0", CELLS[1]), ("C_at_P0_R0", CELLS[2]))
PINS = {
    "counts": ("native_qfo_three_cell_strata_20261007_v1/report.json",
               "55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5"),
    "counts_reader": ("native_qfo_three_cell_strata_readback_20261007_v2.json",
                      "6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969"),
    "protocol": ("SWISS_MODEL_DIVERGENCE_PROTOCOL_20261007.md",
                 "a596371c443471585b9a277cbdd1466bb3f9870200f3187d20286461d3699bfd"),
}
SCOPE = dict(new_bootstrap_draws=0, new_uncertainty=False, independent_confirmation=False,
             new_accuracy_or_resource_admission=False, publication_ready=False,
             scientific_timings_admitted=False, raw_scorer_repeated=False, alignment_repeated=False)
LIMITATIONS = [
    "Development-exposed descriptive strata, not independent confirmation or causal explanation.",
    "Fixed WAG+G4 and retained alignments; model fit and tree-search uncertainty remain unvalidated.",
    "Distance includes paralogy/sampling/composition/domain/alignment effects, not biological time.",
    "Macro-family precision/recall then harmonic F1, not pooled pairs or mean-family F1.",
    "Two conditional contrasts only; candidate-by-reconciliation interaction is not identified.",
    "No subgroup intervals, significance, new bootstrap draws, tuning or default change.",
    "Initial HMM search is on and downstream profile refinement is off in these three cells.",
    "Cell7's failed timing remains ineligible; new construction cost is not inference timing.",
    "Shared-host timing distortion is unknown and potentially tool-dependent.",
]


def need(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve(strict=True)
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def family_statistics(counts):
    need(set(counts) == {"TP", "FP", "FN", "TN"}
         and all(type(v) is int and v >= 0 for v in counts.values())
         and sum(counts.values()) > 0, "Invalid family counts")
    tp, fp, fn = (counts[k] / 2 + 1 for k in ("TP", "FP", "FN"))
    p, r = tp / (tp + fp), tp / (tp + fn)
    return dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)


def aggregate(values):
    if not values:
        return dict.fromkeys(METRICS)
    p, r = (math.fsum(v[m] for v in values) / len(values) for m in ("PPV", "TPR"))
    return dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)


def project(features, feature_reader, counts, count_reader):
    need(features["schema"] == "swiss_model_divergence_features_v1"
         and features["status"] == "features_constructed_unverified"
         and features["failed_families"] == []
         and features["model"] == "WAG+G4" and features["seed"] == 20261007
         and features["prediction_statistics_evaluated"] is False, "Wrong feature scope")
    need(feature_reader["schema"] == "swiss_model_divergence_edge_readback_v1"
         and feature_reader["status"] == "features_verified"
         and feature_reader["families_checked"] == 18 and feature_reader["proteins_checked"] == 563
         and feature_reader["scheduler"]["state"] == "COMPLETED"
         and feature_reader["scheduler"]["exit_code"] == "0:0"
         and features["strata"] == feature_reader["strata"], "Incomplete independent feature verification")
    memberships = features["memberships"]
    families = sorted(memberships)
    genes = [g for f in families for g in memberships[f]]
    need(len(families) == 18 and len(genes) == len(set(genes)) == 563
         and counts["memberships"] == memberships
         and counts["schema"] == "native_qfo_three_cell_strata_v1", "Changed family/member universe")
    need(count_reader["family_rows_checked"] == 54 and count_reader["score_rows_checked"] == 60
         and count_reader["differences_checked"] == 40
         and count_reader["inherited_score_rows_reproduced"] == 40, "Incomplete inherited count reader")
    for doc in (features, feature_reader, counts, count_reader):
        need(all(doc[k] is False for k in ("independent_confirmation", "publication_ready",
                 "new_accuracy_or_resource_admission", "new_uncertainty", "scientific_timings_admitted"))
             and type(doc["new_bootstrap_draws"]) is int and doc["new_bootstrap_draws"] == 0,
             "Inflated inherited scope")
    family_rows = counts["family_rows"]
    need([(r["cell"], r["family"]) for r in family_rows] ==
         [(c, f) for c in CELLS for f in families], "Incomplete retained native counts")
    values = {}
    for row in family_rows:
        value = family_statistics(row["counts_without_prior"])
        need(all(math.isclose(row[m], value[m], rel_tol=0, abs_tol=1e-12) for m in METRICS),
             "Incorrect retained family statistics")
        values[row["cell"], row["family"]] = value
    strata = features["strata"]
    need(list(sorted(strata)) == ["all", "higher_than_median", "lower_or_equal_median"]
         and strata["all"] == families
         and sorted(strata["lower_or_equal_median"] + strata["higher_than_median"]) == families
         and not set(strata["lower_or_equal_median"]) & set(strata["higher_than_median"]),
         "Wrong fixed distance partition")
    rows, differences = [], []
    for name in ("all", "lower_or_equal_median", "higher_than_median"):
        members = strata[name]
        stats = {}
        for cell in CELLS:
            stats[cell] = aggregate([values[cell, f] for f in members])
            rows.append(dict(cell=cell, stratum=name, families=len(members), family_members=members,
                             status="descriptive" if members else "empty",
                             prediction_semantics="native_pair" if cell == CELLS[1] else "group_clique",
                             **stats[cell]))
        for contrast, cell in CONTRASTS:
            delta = {m + "_pp": 100 * (stats[cell][m] - stats[CELLS[0]][m]) if members else None
                     for m in METRICS}
            differences.append(dict(contrast=contrast, candidate=cell, reference=CELLS[0],
                                    stratum=name, families=len(members), family_members=members,
                                    status="descriptive" if members else "empty", **delta))
    for row in rows[:3]:
        original = [r for r in counts["rows"] if r["suite"] == "sequence"
                    and r["stratum"] == "all" and r["cell"] == row["cell"]]
        need(len(original) == 1 and all(math.isclose(row[m], original[0][m], rel_tol=0, abs_tol=1e-12)
                                       for m in METRICS), "All-family statistic changed")
    return dict(schema="native_qfo_swiss_model_divergence_strata_v1", model=features["model"],
                seed=features["seed"], unit=features["unit"], memberships=memberships,
                median_family_distance=features["median_family_distance"], bins=strata,
                features=feature_reader["features"], cells=counts["cells"], family_rows=family_rows,
                rows=rows, differences=differences, limitations=LIMITATIONS, **SCOPE)


def table(report):
    lines = ["# SwissTrees Fixed-Model Divergence Strata", "",
             "Descriptive, development-exposed; WAG+G4 model-estimated distances, not biological time.",
             "Scores are percentages; differences are percentage points. No subgroup confidence intervals.",
             "Macro-family precision/recall then harmonic F1; initial HMM search on, profile refinement off.",
             "", "Median family-distance cutoff: " + str(report["median_family_distance"]), "",
             "| Stratum | Families | Cell | F1 (%) | Precision (%) | Recall (%) |",
             "| --- | ---: | --- | ---: | ---: | ---: |"]
    for row in report["rows"]:
        values = ["NA" if row[m] is None else f"{100 * row[m]:.3f}" for m in METRICS]
        lines.append("| " + " | ".join([row["stratum"], str(row["families"]), row["cell"], *values]) + " |")
    lines.extend(["", "| Stratum | Families | Conditional Contrast | Delta F1 (pp) | Delta Precision (pp) | Delta Recall (pp) |",
                  "| --- | ---: | --- | ---: | ---: | ---: |"])
    for row in report["differences"]:
        values = ["NA" if row[m + "_pp"] is None else f"{row[m + '_pp']:+.3f}" for m in METRICS]
        lines.append("| " + " | ".join([row["stratum"], str(row["families"]), row["contrast"], *values]) + " |")
    lines.extend(["", "## Limitations", "", *("- " + line for line in LIMITATIONS), ""])
    return "\n".join(lines)


def export(repo, features, features_sha, feature_reader, reader_sha, output, source_commit):
    output, repo = Path(output).resolve(), Path(repo).resolve()
    need(not output.exists() and not output.is_symlink(), "Existing projection; never overwrite")
    docs, refs = {}, {}
    for key, path, sha in [
        *((k, repo / "benchmark_tools/results" / name, sha) for k, (name, sha) in PINS.items()),
        ("features", features, features_sha), ("feature_reader", feature_reader, reader_sha),
    ]:
        ref = record(path)
        need(ref["sha256"] == sha, "Changed prospective projection input")
        refs[key] = ref
        if Path(path).suffix == ".json":
            docs[key] = json.loads(Path(path).read_text())
    need(docs["feature_reader"]["report"] == refs["features"]
         and docs["counts_reader"]["report"] == refs["counts"], "Reader/report binding mismatch")
    sources = []
    for doc in docs.values():
        ref = doc["source"]
        need(record(ref["path"]) == ref, "Changed directly inherited source")
        sources.append(ref)
    source = record(__file__)
    need(Path(__file__).read_bytes() == subprocess.check_output(["git", "show",
         source_commit + ":benchmark_tools/export_swiss_model_divergence_strata.py"], cwd=repo),
         "Projection source not committed before execution")
    report = project(docs["features"], docs["feature_reader"], docs["counts"], docs["counts_reader"])
    report.update(source=source, source_commit=source_commit, inputs=refs, inherited_sources=sources)
    output.mkdir()
    outputs = {}
    for key, fields in (("rows", ("stratum", "cell", "families", "status", *METRICS, "prediction_semantics")),
                        ("differences", ("stratum", "contrast", "families", "status", *(m + "_pp" for m in METRICS)))):
        path = output / ("scores.tsv" if key == "rows" else "differences.tsv")
        with path.open("x", encoding="ascii", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(report[key])
        outputs[key] = record(path)
    path = output / "TABLE.md"
    with path.open("x", encoding="ascii") as stream:
        stream.write(table(report))
    outputs["table"] = record(path)
    report["outputs"] = outputs
    with (output / "report.json").open("x", encoding="ascii") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--features-sha256", required=True)
    parser.add_argument("--feature-reader", type=Path, required=True)
    parser.add_argument("--reader-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    args = parser.parse_args()
    result = export(args.repo, args.features, args.features_sha256, args.feature_reader,
                    args.reader_sha256, args.output, args.source_commit)
    print(json.dumps(dict(score_rows=len(result["rows"]), differences=len(result["differences"]))))

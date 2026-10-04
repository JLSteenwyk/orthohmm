"""Portable checked-count replay; never reads historical paths in input metadata."""

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sys


INPUTS = {
    "oracle": ("simulation_gene_tree_oracle_readback_20261004.json", "db016a1b4c68ab39caed061122e2196aadfb3ebdc38a4d4f5ad38d9edd664dec"),
    "upstream": ("simulation_upstream_trace_readback_20261004.json", "729636edf23d1858376860465e0b5f283c6031b301f385cc756cec66f1ed759d"),
    "residual": ("simulation_oracle_residual_readback_20261004.json", "5623a34ff8c7e53aa6f010264d5f9d9521d0023582e9182c8a0b727dda0d607b"),
}
CONDITIONS = {"baseline", "divergent", "turnover", "divergent_turnover",
              "missing20", "uneven_taxa", "taxon_count_control"}
ARMS = ("inferred", "generating_root", "generating_rerooted")
STATUSES = ("oracle_eligible", "mixed_ancestry_ineligible", "unambiguous_bypass")


def require(value, message):
    if not value:
        raise ValueError(message)


def counts(score):
    require(all(type(score[k]) is int and score[k] >= 0 for k in ("tp", "fp", "fn")), "Invalid pair counts")
    tp, fp, fn = [score[k] for k in ("tp", "fp", "fn")]
    values = {"f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0,
              "precision": tp / (tp + fp) if tp + fp else 0.0,
              "recall": tp / (tp + fn) if tp + fn else 0.0}
    for key, value in values.items():
        if key in score:
            require(math.isfinite(score[key]) and math.isclose(score[key], value, rel_tol=0, abs_tol=1e-12),
                    "Metric differs from counts")
    return values


def panel(rows):
    expected = {(c, s) for c in CONDITIONS for s in range(20261101, 20261111)}
    require(len(rows) == 70 and {(r["condition"], r["seed"]) for r in rows} == expected,
            "Incomplete or duplicate 70-cell panel")
    require(all(r["label"] == f"{r['condition']}_{r['seed']}" for r in rows), "Cell label mismatch")
    return {r["label"]: r for r in rows}


def summarize(oracle, upstream, residual):
    require(oracle["status"] == "independent_count_and_candidate_decomposition_verified",
            "Oracle readback is not admitted")
    require(upstream["all_stage_flags_verified_against_native_inputs"] is True,
            "Upstream stage flags were not admitted")
    require(residual["status"] == "independent_xml_rules_and_stage_readback_passed",
            "Residual readback is not admitted")
    o_cells, u_cells = panel(oracle["cells"]), panel(upstream["cells"])
    expected_residuals = {}
    for label, cell in o_cells.items():
        cross = cell["true_pairs_across_candidates"]
        require(type(cross) is int and cross >= 0, "Invalid cross-candidate count")
        for arm in ARMS:
            counts(cell["arms"][arm])
            local = cell["residual_by_arm"][arm]
            require(set(local) == set(STATUSES), "Missing eligibility category")
            for score in local.values():
                counts(score)
            total = {k: sum(s[k] for s in local.values()) for k in ("tp", "fp", "fn")}
            total["fn"] += cross
            require(all(total[k] == cell["arms"][arm][k] for k in total), "Candidate/global decomposition differs")
        expected_residuals[label] = {s: {k: cell["residual_by_arm"]["generating_root"][s][k]
                                        for k in ("fp", "fn")} for s in STATUSES}
        stage = u_cells[label]["totals"]
        require(all(type(v) is int and v >= 0 for v in stage.values()), "Invalid stage count")
        inferred = cell["arms"]["inferred"]
        require(stage.get("native_tp", 0) == inferred["tp"] and stage.get("native_fn", 0) == inferred["fn"]
                and stage.get("true_pairs", 0) == inferred["tp"] + inferred["fn"]
                and stage.get("across_candidates", 0) == cross
                and stage.get("native_fn_within_candidates", 0) + cross == inferred["fn"], "Native/upstream counts differ")
        require(stage.get("different_graph_components", 0) + stage.get("connected_but_separated", 0) == cross,
                "Graph-component partition does not cover cross-candidate losses")
    require(len(oracle["summary"]) == 7 and {r["condition"] for r in oracle["summary"]} == CONDITIONS,
            "Incomplete oracle condition summary")
    require(len(upstream["summary"]) == 7 and {r["condition"] for r in upstream["summary"]} == CONDITIONS,
            "Incomplete upstream condition summary")
    conditions = []
    for condition in sorted(CONDITIONS):
        cells = [r for r in o_cells.values() if r["condition"] == condition]
        means = {arm: {metric: sum(counts(r["arms"][arm])[metric] for r in cells) / len(cells)
                       for metric in ("f1", "precision", "recall")} for arm in ARMS}
        declared = next(r for r in oracle["summary"] if r["condition"] == condition)
        require(declared["scored_cells"] == 10 and all(math.isclose(means[a][m], declared["mean_metrics"][a][m],
                rel_tol=0, abs_tol=1e-12) for a in ARMS for m in means[a]), "Oracle condition mean mismatch")
        cross = sum(r["true_pairs_across_candidates"] for r in cells)
        fn = sum(r["arms"]["generating_root"]["fn"] for r in cells)
        stage = dict(sum((Counter(u_cells[r["label"]]["totals"]) for r in cells), Counter()))
        u_declared = next(r for r in upstream["summary"] if r["condition"] == condition)
        require(stage == u_declared["totals"] and u_declared["cells"] == 10, "Upstream summary projection mismatch")
        require(fn == declared["oracle_total_fn"] and cross == declared["true_pairs_across_candidates"],
                "Oracle condition error projection mismatch")
        conditions.append({"condition": condition, "cells": 10, "mean_metrics": means,
            "generating_root_f1_change_pp": 100 * (means["generating_root"]["f1"] - means["inferred"]["f1"]),
            "oracle_fn": fn, "cross_candidate_fn": cross,
            "cross_candidate_oracle_fn_share": cross / fn if fn else None,
            "upstream_totals": stage})
    seen, status_counts, error_classes = set(), Counter(), Counter()
    observed_residuals = {label: {s: {"fp": 0, "fn": 0} for s in STATUSES} for label in o_cells}
    for row in residual["candidates"]:
        key = (row["cell"], row["family"])
        require(key not in seen and row["cell"] in o_cells and row["status"] in STATUSES, "Invalid residual cohort")
        seen.add(key)
        counts(row["counts"])
        require(row["counts"]["fp"] or row["counts"]["fn"], "Residual candidate has no error")
        require(sum(row["error_classes"].values()) == row["counts"]["fp"] + row["counts"]["fn"],
                "Residual classes do not cover errors")
        for k in ("fp", "fn"):
            observed_residuals[row["cell"]][row["status"]][k] += row["counts"][k]
        status_counts[row["status"]] += 1
        error_classes.update(row["error_classes"])
    require(observed_residuals == expected_residuals, "Residual cohort does not account for all original within-candidate errors")
    summary = {"candidates": len(seen), "status_counts": dict(status_counts),
        "pair_rows": sum(r["pair_rows_verified"] for r in residual["candidates"]),
        "counts": {k: sum(r["counts"][k] for r in residual["candidates"]) for k in ("tp", "fp", "fn")},
        "error_classes": dict(error_classes)}
    require(summary == residual["summary"], "Residual summary differs from candidates")
    require(sum(r["totals"]["true_pairs"] for r in u_cells.values()) == upstream["pair_rows_verified"],
            "Upstream all-pair count mismatch")
    return {"conditions": conditions, "screened_cells": 70,
        "screened_candidates": sum(sum(r["candidate_counts"].values()) for r in o_cells.values()),
        "eligible_candidates": sum(r["candidate_counts"].get("oracle_eligible", 0) for r in o_cells.values()),
        "topology": oracle["topology"], "upstream_true_pair_rows": upstream["pair_rows_verified"],
        "residual": summary}


def render(value):
    lines = ["# Simulation Mechanism Summary", "",
        "Finite-panel reporting replay from checked counts; not raw admission, native inference, confidence intervals or independent validation.", "",
        "| Condition | Inferred F1 (%) | Generating-root F1 (%) | Change (pp) | Oracle FN | Cross-candidate FN | Different graph components | Connected but separated |",
        "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for row in value["conditions"]:
        stage, means = row["upstream_totals"], row["mean_metrics"]
        lines.append(f"| {row['condition']} | {100*means['inferred']['f1']:.3f} | {100*means['generating_root']['f1']:.3f} | "
            f"{row['generating_root_f1_change_pp']:+.3f} | {row['oracle_fn']:,} | {row['cross_candidate_fn']:,} | "
            f"{stage.get('different_graph_components', 0):,} | {stage.get('connected_but_separated', 0):,} |")
    lines.extend(["", "## Complete Within-Candidate Residual Cohort", "", "| Mechanism | Error Pairs |", "|---|---:|"])
    for key, count in sorted(value["residual"]["error_classes"].items()):
        lines.append(f"| {key} | {count} |")
    lines.extend(["", f"Screened {value['screened_cells']} cells / {value['screened_candidates']:,} candidates; "
        f"{value['eligible_candidates']:,} eligible tree controls. Upstream: {value['upstream_true_pair_rows']:,} true-pair rows.",
        f"Residual cohort: {value['residual']['candidates']} candidates / {value['residual']['pair_rows']} pair rows.",
        "Missing significant hits do not isolate search-stage causes; candidate paths may cross other families.",
        "Retained candidate membership and constraint policy are fixed. No defaults or original scores change.", ""])
    return "\n".join(lines)


def run(directory, output):
    require(not output.exists() and not output.is_symlink(), "Refusing existing output")
    reports = {}
    for role, (name, digest) in INPUTS.items():
        data = (directory / name).read_bytes()
        require(hashlib.sha256(data).hexdigest() == digest, "Changed checked input: " + name)
        reports[role] = json.loads(data)
    result = {"status": "checked_simulation_mechanism_reporting_replayed",
              "summary": summarize(**reports), "input_sha256": {name: digest for name, digest in INPUTS.values()},
              "python": sys.version, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "historical_native_paths_accessed": False, "native_inference_repeated": False,
              "publication_ready": False, "scope": "Derived checked-summary arithmetic only"}
    output.mkdir(parents=True)
    (output / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    (output / "summary.md").write_text(render(result["summary"]))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.inputs.resolve(), args.output.absolute())
    print(json.dumps({"status": result["status"], "cells": result["summary"]["screened_cells"],
                      "residual": result["summary"]["residual"]}, sort_keys=True))

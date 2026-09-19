"""Independently reconstruct frozen parameter contrasts from admitted counts."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

ARMS = ("control", "cpm_low", "cpm_high", "norm_low", "norm_high", "margin_low", "margin_high")
METRICS = ("F1", "PPV", "TPR")


def verify(report):
    if (report["status"] != "corrected_qfo_parameter_uncertainty_audited"
            or report["scientific_inputs_admitted"] is not True
            or report["publication_ready"] is not False
            or report["replicates"] != 100000 or report["seed"] != 20260925
            or report["multiplicity_endpoints"] != 18 or report["alpha"] != .05
            or report["quantile_method"] != "linear" or report["protocol_controls_match"] is not True):
        raise ValueError("Wrong result status or frozen controls")
    counts = report["reconstructed_counts"]
    families = report["families"]
    rows = counts["arms"]
    if (counts["status"] != "corrected_qfo_parameter_swiss_counts_verified"
            or counts["publication_ready"] is not False or counts["uncertainty_admitted"] is not False
            or counts["shared_represented_genes"] or counts["families"] != families
            or len(families) != 18 or len(set(families)) != 18
            or tuple(r["arm"] for r in rows) != ARMS
            or set(report["point_estimates"]) != set(ARMS)):
        raise ValueError("Wrong admitted arm or family inventory")
    if report["arms"] != [{k: r[k] for k in ("arm", "status", "reason") if k in r} for r in rows]:
        raise ValueError("Arm availability summary differs")

    def scores(pr):
        p, r = pr[..., 0], pr[..., 1]
        return np.stack((2 * p * r / (p + r), p, r), axis=-1)

    def compare(actual, expected):
        actual, expected = np.asarray(actual), np.asarray(expected)
        if actual.shape != expected.shape or not np.allclose(actual, expected, rtol=0, atol=1e-12):
            raise ValueError("Parameter arithmetic differs from independent reproduction")

    weights = np.random.Generator(np.random.PCG64(20260925)).multinomial(18, [1 / 18] * 18, size=100000)
    points, draws, family_scores, universe = {}, {}, {}, None
    for row in rows:
        name = row["arm"]
        if row["status"] == "not_admitted":
            if (not isinstance(row.get("reason"), str) or not row["reason"].strip()
                    or "families" in row or "aggregate" in row or report["point_estimates"][name] is not None):
                raise ValueError("Imputed or unexplained unavailable arm")
            continue
        if row["status"] != "counts_verified" or [r["family"] for r in row["families"]] != families:
            raise ValueError("Wrong admitted family inventory")
        pr, truth, seen, total = [], [], set(), 0
        for family in row["families"]:
            raw, genes = family["counts_without_prior"], family["represented_genes"]
            if (set(raw) != {"TP", "FP", "FN", "TN"}
                    or any(type(v) is not int or v < 0 for v in raw.values())
                    or not sum(raw.values()) or len(genes) <= 5
                    or len(set(genes)) != len(genes) or seen.intersection(genes)):
                raise ValueError("Invalid counts or overlapping family members")
            seen.update(genes)
            truth.append((sorted(genes), raw["TP"] + raw["FN"], raw["FP"] + raw["TN"]))
            total += sum(raw.values())
            tp, fp, fn = (raw[k] / 2 + 1 for k in ("TP", "FP", "FN"))
            pr.append((tp / (tp + fp), tp / (tp + fn)))
            compare([family["statistics_with_prior"][k] for k in METRICS], scores(np.asarray(pr[-1])))
        if universe is None:
            universe = truth
        if truth != universe or total != counts["reference_relation_count"]:
            raise ValueError("Reference truth or membership differs across arms")
        pr = np.asarray(pr)
        points[name], family_scores[name] = scores(pr.mean(axis=0)), scores(pr)
        compare([row["aggregate"][k] for k in METRICS], points[name])
        compare([report["point_estimates"][name][k] for k in METRICS], points[name])
        means = np.zeros((100000, 2))
        for i in range(18):
            means += pr[i] * weights[:, i, None] / 18
        draws[name] = scores(means)

    if len(report["comparisons"]) != 6:
        raise ValueError("Wrong contrast count")
    estimated = 0
    for row, candidate in zip(report["comparisons"], ARMS[1:]):
        if (row["candidate"], row["reference"]) != (candidate, "control"):
            raise ValueError("Wrong contrast identity")
        if candidate not in points or "control" not in points:
            reason = "baseline_not_admitted" if "control" not in points else "variant_not_admitted"
            if (row["status"] != "not_estimable" or row["reason"] != reason
                    or row["metrics"] is not None or row["family_differences"] is not None):
                raise ValueError("Unavailable contrast was not preserved")
            continue
        if row["status"] != "estimated":
            raise ValueError("Available contrast omitted")
        estimated += 1
        sample = draws[candidate] - draws["control"]
        delta = family_scores[candidate] - family_scores["control"]
        if [r["family"] for r in row["family_differences"]] != families:
            raise ValueError("Wrong per-family contrast inventory")
        compare([[r[m] for m in METRICS] for r in row["family_differences"]], delta)
        for j, metric in enumerate(METRICS):
            item = row["metrics"][metric]
            compare(item["difference"], points[candidate][j] - points["control"][j])
            for key, q in (("paired_percentile_ci", [.025, .975]),
                           ("bonferroni_percentile_ci", [.05 / 36, 1 - .05 / 36])):
                compare(item[key], np.quantile(sample[:, j], q, method="linear"))
            expected = [int(np.sum(delta[:, j] > 1e-10)), int(np.sum(np.abs(delta[:, j]) <= 1e-10)),
                        int(np.sum(delta[:, j] < -1e-10))]
            if [item[k] for k in ("family_wins", "family_ties", "family_losses")] != expected:
                raise ValueError("Wrong family wins/ties/losses")
    if (report["estimated_contrasts"] != estimated or report["complete_panel"] is not (estimated == 6)
            or report["uncertainty_admitted"] is not (estimated > 0)):
        raise ValueError("Incorrect uncertainty admission or completeness")
    return estimated * 3


def run(path, digest, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    result = read_frozen(path, digest)
    source = record(__file__)
    checked = [source, record(path), result["source"], *result["helpers"], *result["checked_inputs"]]
    for item in checked:
        check(item)
    endpoints = verify(result)
    for item in checked:
        check(item)
    report = {"status": "qfo_parameter_uncertainty_numerically_reproduced", "endpoints": endpoints,
        "planned_endpoints": 18, "input": record(path), "source": source, "numpy_version": np.__version__,
        "absolute_tolerance": 1e-12, "publication_ready": False,
        "limitations": ["Reconstructs arithmetic from admitted counts, not raw predictions or native inference.",
                        "Uses the same NumPy random generator and quantile implementation.",
                        "Zero endpoints means a control-only or unavailable-control check, not parameter robustness."]}
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--results-sha256", dest="digest", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.results.resolve(), args.digest, args.output.absolute())

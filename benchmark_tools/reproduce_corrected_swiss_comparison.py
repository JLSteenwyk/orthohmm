"""Independently reproduce corrected comparator arithmetic from admitted counts."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

METHODS = ("orthohmm_high_sensitivity", "orthohmm_phylogeny_satellite_v2", "orthofinder_3_1_5_full",
           "orthofinder_3_1_5_sequence_only", "sonicparanoid_2_0_9", "proteinortho_6_3_6",
           "fastoma_0_3_5", "orthomcl_1_4")
CONTRASTS = ((0, 2), (1, 2), (3, 2), (4, 2), (5, 2), (6, 2), (7, 2), (1, 0))
METRICS = ("F1", "PPV", "TPR")


def verify(report):
    if (report["status"] != "corrected_swiss_comparison_intervals_audited"
            or report["scientific_inputs_admitted"] is not True or report["publication_ready"] is not False
            or report["replicates"] != 100000 or report["seed"] != 20260920
            or report["multiplicity_endpoints"] != 24 or report["alpha"] != .05
            or report["quantile_method"] != "linear" or report["protocol_controls_match"] is not True):
        raise ValueError("Wrong result status or frozen controls")
    rows = report["reconstructed_counts"]["methods"]
    families = report["families"]
    if tuple(r["method"] for r in rows) != METHODS or len(families) != 18 or len(set(families)) != 18:
        raise ValueError("Wrong method or family inventory")
    weights = np.random.Generator(np.random.PCG64(20260920)).multinomial(18, [1 / 18] * 18, size=100000)
    points, draws, family_scores = {}, {}, {}

    def scores(pr):
        p, r = pr[..., 0], pr[..., 1]
        return np.stack((2 * p * r / (p + r), p, r), axis=-1)

    def compare(actual, expected):
        if not np.allclose(actual, expected, rtol=0, atol=1e-12):
            raise ValueError("Comparator arithmetic differs from independent reproduction")

    for row in rows:
        name = row["method"]
        if row["status"] == "not_admitted":
            if report["point_estimates"][name] is not None:
                raise ValueError("Imputed unavailable method")
            continue
        if row["status"] != "counts_verified" or [r["family"] for r in row["families"]] != families:
            raise ValueError("Wrong admitted family inventory")
        pr = []
        for family in row["families"]:
            raw = family["counts_without_prior"]
            tp, fp, fn = (raw[k] / 2 + 1 for k in ("TP", "FP", "FN"))
            pr.append((tp / (tp + fp), tp / (tp + fn)))
        pr = np.asarray(pr)
        points[name], family_scores[name] = scores(pr.mean(axis=0)), scores(pr)
        compare([report["point_estimates"][name][k] for k in METRICS], points[name])
        means = np.zeros((100000, 2))
        for i in range(18):
            means += pr[i] * weights[:, i, None] / 18
        draws[name] = scores(means)
    if len(report["comparisons"]) != 8:
        raise ValueError("Wrong contrast inventory")
    estimated = 0
    for row, (a, b) in zip(report["comparisons"], CONTRASTS):
        candidate, reference = METHODS[a], METHODS[b]
        if (row["candidate"], row["reference"]) != (candidate, reference):
            raise ValueError("Wrong contrast identity")
        missing = [name for name in (candidate, reference) if name not in points]
        if missing:
            if (row["status"] != "not_estimable" or row["unavailable_methods"] != missing
                    or row["metrics"] is not None or row["family_differences"] is not None):
                raise ValueError("Unavailable contrast was not preserved")
            continue
        if row["status"] != "estimated":
            raise ValueError("Available contrast omitted")
        estimated += 1
        point = points[candidate] - points[reference]
        sample = draws[candidate] - draws[reference]
        deltas = family_scores[candidate] - family_scores[reference]
        if [r["family"] for r in row["family_differences"]] != families:
            raise ValueError("Wrong per-family contrast inventory")
        compare([[r[m] for m in METRICS] for r in row["family_differences"]], deltas)
        for j, metric in enumerate(METRICS):
            item = row["metrics"][metric]
            compare(item["difference"], point[j])
            for key, quantiles in (("paired_percentile_ci", [.025, .975]),
                                   ("bonferroni_percentile_ci", [.05 / 48, 1 - .05 / 48])):
                compare(item[key], np.quantile(sample[:, j], quantiles, method="linear"))
            expected = [int(np.sum(deltas[:, j] > 1e-10)), int(np.sum(np.abs(deltas[:, j]) <= 1e-10)),
                        int(np.sum(deltas[:, j] < -1e-10))]
            if [item[k] for k in ("family_wins", "family_ties", "family_losses")] != expected:
                raise ValueError("Wrong family wins/ties/losses")
    if (report["estimated_contrasts"] != estimated or report["complete_panel"] is not (estimated == 8)
            or report["uncertainty_admitted"] is not (estimated > 0)):
        raise ValueError("Incorrect interval admission or completeness")
    return estimated * 3


def run(path, digest, output):
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.run_simulation_methods import read_frozen

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
    report = {"status": "corrected_swiss_comparison_numerically_reproduced", "endpoints": endpoints,
        "planned_endpoints": 24, "input": record(path), "source": source, "numpy_version": np.__version__,
        "absolute_tolerance": 1e-12, "publication_ready": False,
        "limitations": ["Recomputes arithmetic from admitted counts, not an independent raw-count or inference audit.",
                        "Uses the same NumPy random generator and quantile implementation."]}
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


def retained_counts_only(path, digest, output):
    """Portable arithmetic check; never follow historical source/evidence paths."""
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != digest:
        raise ValueError("Retained counts report hash differs")
    source = Path(__file__).read_bytes()
    endpoints = verify(json.loads(content))
    if path.read_bytes() != content or Path(__file__).read_bytes() != source:
        raise ValueError("Input or verifier changed during reproduction")
    report = dict(status="retained_corrected_swiss_arithmetic_reproduced",
        endpoints=endpoints, planned_endpoints=24, input_sha256=digest,
        source_sha256=hashlib.sha256(source).hexdigest(), python_version=sys.version,
        numpy_version=np.__version__, absolute_tolerance=1e-12,
        raw_source_admission_repeated=False, historical_paths_accessed=False,
        publication_ready=False,
        limitations=["Arithmetic only from hash-pinned retained counts; native inference, conversion and raw scoring were not rerun.",
            "Historical absolute paths are provenance, not accessed inputs.",
            "Uses the same NumPy RNG and quantile implementation; this is not cross-platform validation or a full release."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--results-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--retained-counts-only", action="store_true",
                        help="Recheck portable arithmetic without repeating historical source admission")
    args = parser.parse_args()
    runner = retained_counts_only if args.retained_counts_only else run
    runner(args.results.resolve(), args.results_sha256, args.output.absolute())

"""Check strata arithmetic independently of the production bootstrap helpers."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

METHODS = ("high_sensitivity", "phylogenetic", "orthofinder_full")
METRICS = ("F1", "PPV", "TPR")
CONTRASTS = ((0, 2), (1, 2), (1, 0))


def verify(report):
    if (report["status"] != "corrected_swiss_primary_stratified_intervals"
            or report["scientific_inputs_admitted"] is not True or report["uncertainty_admitted"] is not True
            or report["replicates"] != 100000 or report["seed"] != 20260924
            or report["multiplicity_endpoints"] != 27 or report["alpha"] != .05
            or report["quantile_method"] != "linear" or report["publication_ready"] is not False):
        raise ValueError("Wrong admitted strata result or frozen numerical settings")
    cells = report["reconstructed_counts"]["factorial"]["cells"]
    if [c["cell"] for c in cells] != [f"p{p}_c{c}_r{r}" for p in (0, 1) for c in (0, 1) for r in (0, 1)]:
        raise ValueError("Wrong factorial cell inventory")
    of = report["reconstructed_counts"]["orthofinder"]
    if of["method"] != "orthofinder_full":
        raise ValueError("Wrong comparator")
    rows = [{r["family"]: r for r in source["families"]} for source in (cells[4], cells[7], of)]
    families = [f for b in ("lower", "higher", "missing") for f in report["bins"][b]["families"]]
    if len(families) != 18 or len(set(families)) != 18 or any(set(r) != set(families) for r in rows):
        raise ValueError("Wrong family inventory")
    rng = np.random.Generator(np.random.PCG64(20260924))
    differences, points, endpoints = {}, {}, 0

    def scores(pr):
        p, r = pr[..., 0], pr[..., 1]
        return np.stack((2 * p * r / (p + r), p, r), axis=-1)

    def compare(actual, expected):
        if not np.allclose(actual, expected, rtol=0, atol=1e-12):
            raise ValueError("Strata arithmetic differs from independent reproduction")

    def interval(row, point, draws):
        if point is None:
            if row["difference"] is not None:
                raise ValueError("Unexpected missing-bin point")
        else:
            compare(row["difference"], point)
        for key, q in (("paired_percentile_ci", [.025, .975]),
                       ("bonferroni_percentile_ci", [.05 / 54, 1 - .05 / 54])):
            if draws is None:
                if row[key] is not None:
                    raise ValueError("Unexpected ineligible interval")
            else:
                compare(row[key], np.quantile(draws, q, method="linear"))

    for name in ("lower", "higher", "missing"):
        bin_result = report["bins"][name]
        selected = bin_result["families"]
        n = len(selected)
        eligible = name != "missing" and n >= 5
        if bin_result["interval_eligible"] is not eligible:
            raise ValueError("Wrong interval eligibility")
        pr = np.empty((3, n, 2))
        for m, method in enumerate(rows):
            for i, family in enumerate(selected):
                raw = method[family]["counts_without_prior"]
                tp, fp, fn = (raw[k] / 2 + 1 for k in ("TP", "FP", "FN"))
                pr[m, i] = (tp / (tp + fp), tp / (tp + fn))
        point = scores(pr.mean(axis=1)) if n else None
        sampled = None
        if eligible:
            weights = rng.multinomial(n, [1 / n] * n, size=100000)
            # Sum each family's weighted contribution without the production matrix product.
            means = np.zeros((3, 100000, 2))
            for i in range(n):
                means += pr[:, i, None, :] * weights[None, :, i, None] / n
            sampled = scores(means)
        points[name], differences[name] = [], []
        for m, method in enumerate(METHODS):
            actual = bin_result["point_estimates"][method]
            if point is None:
                if actual is not None:
                    raise ValueError("Unexpected empty-bin estimate")
            else:
                compare([actual[k] for k in METRICS], point[m])
        if len(bin_result["comparisons"]) != 3:
            raise ValueError("Wrong contrast inventory")
        for row, (a, b) in zip(bin_result["comparisons"], CONTRASTS):
            if (row["candidate"], row["reference"]) != (METHODS[a], METHODS[b]):
                raise ValueError("Wrong contrast identity")
            delta = None if point is None else point[a] - point[b]
            draws = None if sampled is None else sampled[a] - sampled[b]
            points[name].append(delta)
            differences[name].append(draws)
            family_delta = scores(pr[a]) - scores(pr[b])
            for j, metric in enumerate(METRICS):
                item = row["metrics"][metric]
                interval(item, None if delta is None else delta[j], None if draws is None else draws[:, j])
                expected = [int(np.sum(family_delta[:, j] > 1e-10)),
                            int(np.sum(np.abs(family_delta[:, j]) <= 1e-10)),
                            int(np.sum(family_delta[:, j] < -1e-10))]
                if [item[k] for k in ("family_wins", "family_ties", "family_losses")] != expected:
                    raise ValueError("Wrong family win/tie/loss counts")
                endpoints += int(name != "missing")
    if len(report["interactions"]) != 3:
        raise ValueError("Wrong interaction inventory")
    for i, (row, (a, b)) in enumerate(zip(report["interactions"], CONTRASTS)):
        if (row["candidate"], row["reference"]) != (METHODS[a], METHODS[b]):
            raise ValueError("Wrong interaction identity")
        lo, hi = points["lower"][i], points["higher"][i]
        ld, hd = differences["lower"][i], differences["higher"][i]
        for j, metric in enumerate(METRICS):
            interval(row["metrics"][metric], None if lo is None or hi is None else hi[j] - lo[j],
                     None if ld is None or hd is None else hd[:, j] - ld[:, j])
            endpoints += 1
    return endpoints


def run(path, digest, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    result = read_frozen(path, digest)
    inputs = [record(path), result["source"], *result["helpers"], *result["checked_inputs"]]
    for item in inputs:
        check(item)
    endpoints = verify(result)
    for item in inputs:
        check(item)
    report = {"status": "corrected_swiss_strata_numerically_reproduced", "endpoints": endpoints,
        "input": record(path), "source": record(__file__), "numpy_version": np.__version__,
        "absolute_tolerance": 1e-12, "publication_ready": False,
        "limitations": ["Recomputes arithmetic from admitted counts; not a new raw-count or inference audit.",
                        "Uses the same NumPy random generator and quantile implementation."]}
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--results-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.results.resolve(), args.results_sha256, args.output.absolute())

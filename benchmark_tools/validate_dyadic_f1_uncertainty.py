"""Synthetic validation of a candidate dyadic F1 variance; no benchmark CIs."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

SEED = 20260927
SIZES = (64, 256)
CASES = ("independent_sparse", "node_sparse", "node_dense", "clade_sparse")


def f1_gradient(counts):
    counts = np.asarray(counts, dtype=float)
    if counts.shape != (3,) or not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError("Require three finite nonnegative TP/FP/FN counts")
    tp, fp, fn = counts
    denominator = 2 * tp + fp + fn
    if denominator <= 0:
        raise ValueError("Undefined F1")
    return 2 * tp / denominator, np.array([2 * (fp + fn), -2 * tp, -2 * tp]) / denominator**2


def shared_endpoint_variance(single, edge, left, right):
    """Exact overlap-sum algebra for centered linearized contributions.

    Each edge's squared term appears at both endpoints, so subtract it once.
    Diagonal (single-family) terms enter only their own endpoint. Statistical
    validity additionally requires assumptions tested by the simulation.
    """
    total = np.asarray(single, dtype=float).copy()
    np.add.at(total, left, edge)
    np.add.at(total, right, edge)
    return float(total @ total - edge @ edge)


def contrast(diagonal, edges, left, right):
    n = len(diagonal)
    if (diagonal.shape != (n, 2, 3) or edges.shape != (len(left), 2, 3)
            or not np.isfinite(diagonal).all() or not np.isfinite(edges).all()
            or (diagonal < 0).any() or (edges < 0).any()):
        raise ValueError("Invalid paired count arrays")
    expected_left, expected_right = np.triu_indices(n, 1)
    if not (np.array_equal(left, expected_left) and np.array_equal(right, expected_right)):
        raise ValueError("Require the complete unordered family-pair grid, including zeros")
    counts = diagonal.sum(axis=0) + edges.sum(axis=0)
    scores, gradients = zip(*(f1_gradient(c) for c in counts))
    signed = np.array([gradients[0], -gradients[1]])
    d = np.einsum("nmc,mc->n", diagonal, signed)
    e = np.einsum("nmc,mc->n", edges, signed)
    # Diagonal and off-diagonal opportunities have different distributions.
    # Center each stratum separately; including zero dyads is essential.
    d -= d.mean()
    e -= e.mean()
    variance = shared_endpoint_variance(d, e, left, right)
    naive = float(d @ d + e @ e)
    return float(scores[0] - scores[1]), variance, naive


def sample(rng, n, case, left, right):
    if case not in CASES:
        raise ValueError("Unknown simulation case")
    if case == "independent_sparse":
        x = np.zeros(n)
    elif case == "clade_sparse":
        # Adjacent groups of eight share a latent effect: disjoint dyads can
        # now be dependent, deliberately violating the candidate's assumption.
        x = np.repeat(rng.choice([-1.0, 1.0], size=n // 8), 8)
    else:
        x = rng.choice([-1.0, 1.0], size=n)
    probabilities = np.column_stack((0.75 + 0.15 * x, 0.80 + 0.10 * x))
    uniforms = rng.random((n, 20, 1))
    tp = (uniforms < probabilities[:, None, :]).sum(axis=1)
    diagonal = np.zeros((n, 2, 3))
    diagonal[:, :, 0] = tp
    diagonal[:, :, 2] = 20 - tp
    divisor = 32 if case == "node_dense" else n
    rate = (1 + 0.35 * (x[left] + x[right])) / divisor
    common = rng.poisson(5 * rate)
    extra = rng.poisson(3 * rate)
    edges = np.zeros((len(left), 2, 3))
    edges[:, 0, 1] = common + extra
    edges[:, 1, 1] = common
    expected = np.array([[15 * n, 8 * len(left) / divisor, 5 * n],
                         [16 * n, 5 * len(left) / divisor, 4 * n]])
    truth = f1_gradient(expected[0])[0] - f1_gradient(expected[1])[0]
    return diagonal, edges, truth


def wilson(successes, trials):
    p = successes / trials
    z = 1.959963984540054
    denominator = 1 + z * z / trials
    center = (p + z * z / (2 * trials)) / denominator
    half = z * np.sqrt(p * (1 - p) / trials + z * z / (4 * trials**2)) / denominator
    return [float(center - half), float(center + half)]


def run(replicates):
    if replicates < 2:
        raise ValueError("Need at least two replicates")
    rows = []
    for case_id, case in enumerate(CASES):
        for n in SIZES:
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED, case_id, n])))
            left, right = np.triu_indices(n, 1)
            estimates, variances, naive_variances = [], [], []
            covered = naive_covered = failed = 0
            for _ in range(replicates):
                diagonal, edges, truth = sample(rng, n, case, left, right)
                estimate, variance, naive = contrast(diagonal, edges, left, right)
                estimates.append(estimate)
                variances.append(variance)
                naive_variances.append(naive)
                if variance <= 0 or not np.isfinite(variance):
                    failed += 1
                else:
                    covered += abs(estimate - truth) <= 1.959963984540054 * np.sqrt(variance)
                naive_covered += abs(estimate - truth) <= 1.959963984540054 * np.sqrt(naive)
            rows.append(dict(case=case, families=n, replicates=replicates,
                known_ratio_of_expected_counts_difference=float(truth),
                mean_estimate=float(np.mean(estimates)), bias=float(np.mean(estimates) - truth),
                empirical_variance=float(np.var(estimates, ddof=1)),
                mean_candidate_variance=float(np.mean(variances)),
                mean_naive_variance=float(np.mean(naive_variances)),
                invalid_candidate_variances=failed, candidate_covered=int(covered),
                candidate_coverage=float(covered / replicates),
                candidate_coverage_wilson95=wilson(covered, replicates),
                naive_coverage=float(naive_covered / replicates),
                prespecified_screen_pass=bool(failed == 0 and covered / replicates >= 0.925)
                    if case != "clade_sparse" else None))
    return dict(status="synthetic_variance_screen_complete", seed=SEED,
        generator="PCG64 with SeedSequence([seed, case_index, family_count])", numpy_version=np.__version__,
        rows=rows, benchmark_intervals_admitted=False,
        all_in_model_cells_pass=all(r["prespecified_screen_pass"] for r in rows if r["case"] != "clade_sparse"),
        limitations=["Candidate delta-method variance, not a theorem for native VGNC F1.",
            "Synthetic exchangeable families do not establish biological family independence.",
            "Known target is the difference of ratios of expected counts, not the expectation of sample F1.",
            "Invalid variances count as noncoverage; no clipping, retries or case omission.",
            "Coverage screen is necessary exploratory evidence, not authorization to publish benchmark CIs."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = run(args.replicates)
    result["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")

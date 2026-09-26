"""Reference-size and rare-error stress screen; no native benchmark intervals."""

import argparse
import json
from pathlib import Path

import numpy as np

from benchmark_tools.acquire_publication_fasttree import identity, verify
from benchmark_tools.validate_dyadic_f1_uncertainty import f1_gradient, wilson

REFERENCE = "benchmark_tools/results/corrected_vgnc_blocks_20260926.json"
REFERENCE_SHA = "01a51a62a536bcf759f2a0ef83c56196a92022b0bb8290448eacbfefa47d6005"
SIZES = (256, 16844)
CASES = ("unequal_regular", "unequal_node", "rare_perfect")
SEED = 20260928


def sparse_contrast(diagonal, edge_counts, left, right):
    n = len(diagonal)
    if n < 2 or diagonal.shape != (n, 2, 3) or edge_counts.shape != (len(left), 2, 3):
        raise ValueError("Invalid count shapes")
    if (len(right) != len(left) or np.any(left < 0) or np.any(left >= right) or np.any(right >= n)
            or len(np.unique(left * n + right)) != len(left)):
        raise ValueError("Require unique unordered in-range dyads")
    if any(not np.isfinite(v).all() or (v < 0).any() for v in (diagonal, edge_counts)):
        raise ValueError("Invalid counts")
    totals = diagonal.sum(axis=0) + edge_counts.sum(axis=0)
    a, ga = f1_gradient(totals[0])
    b, gb = f1_gradient(totals[1])
    gradient = np.array([ga, -gb])
    d = np.einsum("nmc,mc->n", diagonal, gradient)
    e = np.einsum("nmc,mc->n", edge_counts, gradient)
    d -= d.mean()
    dyads = n * (n - 1) // 2
    mean = e.sum() / dyads
    # Include omitted zero dyads analytically, not as absent opportunities.
    incident = d - (n - 1) * mean
    np.add.at(incident, left, e)
    np.add.at(incident, right, e)
    edge_squares = e @ e - dyads * mean**2
    return float(a - b), float(incident @ incident - edge_squares)


def poisson_edges(rng, n, total_mean, x):
    # Uniform proposals over ordered distinct endpoints give uniform unordered
    # dyads. Thinning yields the specified additive endpoint-effect intensity.
    bound = 1.7 if np.any(x) else 1.
    count = rng.poisson(total_mean * bound)
    left = rng.integers(n, size=count)
    right = rng.integers(n - 1, size=count)
    right += right >= left
    keep = rng.random(count) < (1 + .35 * (x[left] + x[right])) / bound
    left, right = np.minimum(left[keep], right[keep]), np.maximum(left[keep], right[keep])
    return left * n + right


def sample(rng, n, case, sizes, probabilities):
    if case not in CASES:
        raise ValueError("Unknown case")
    trials = rng.choice(sizes, size=n, p=probabilities)
    x = rng.choice([-1., 1.], size=n) if case == "unequal_node" else np.zeros(n)
    recall_a, recall_b = (.75 + .15 * x, .80 + .10 * x)
    means = (4. * (n - 1), 2.5 * (n - 1))
    if case == "rare_perfect":
        recall_a, recall_b = np.ones(n), np.ones(n)
        means = (2., 1.)
    tp_b = rng.binomial(trials, recall_b)
    tp_a = rng.binomial(tp_b, recall_a / recall_b)
    diagonal = np.zeros((n, 2, 3))
    diagonal[:, :, 0] = np.column_stack((tp_a, tp_b))
    diagonal[:, :, 2] = trials[:, None] - diagonal[:, :, 0]
    common = poisson_edges(rng, n, means[1], x)
    extra = poisson_edges(rng, n, means[0] - means[1], x)
    keys, inverse = np.unique(np.concatenate((common, extra)), return_inverse=True)
    edge_counts = np.zeros((len(keys), 2, 3))
    edge_counts[:, 0, 1] = np.bincount(inverse, minlength=len(keys))
    edge_counts[:, 1, 1] = np.bincount(inverse[:len(common)], minlength=len(keys))
    expected_trials = n * float(sizes @ probabilities)
    recalls = (1., 1.) if case == "rare_perfect" else (.75, .80)
    targets = [f1_gradient([expected_trials * r, f, expected_trials * (1 - r)])[0]
               for r, f in zip(recalls, means)]
    return diagonal, edge_counts, keys // n, keys % n, targets[0] - targets[1]


def run(repo, replicates):
    if replicates < 2:
        raise ValueError("Need at least two replicates")
    reference = verify(repo / REFERENCE, REFERENCE_SHA)
    raw = json.loads((repo / REFERENCE).read_text())
    histogram = raw["reference_pair_count_histogram"]
    sizes = np.array(sorted(map(int, histogram)))
    weights = np.array([histogram[str(k)] for k in sizes])
    if weights.sum() != 16844 or sizes @ weights != 23934:
        raise ValueError("Changed reference-only size distribution")
    probabilities = weights / weights.sum()
    rows = []
    for case_id, case in enumerate(CASES):
        for n in SIZES:
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED, case_id, n])))
            estimates, variances = [], []
            covered = invalid = 0
            for _ in range(replicates):
                diagonal, edges, left, right, truth = sample(rng, n, case, sizes, probabilities)
                estimate, variance = sparse_contrast(diagonal, edges, left, right)
                estimates.append(estimate)
                variances.append(variance)
                if variance <= 0 or not np.isfinite(variance):
                    invalid += 1
                else:
                    covered += abs(estimate - truth) <= 1.959963984540054 * np.sqrt(variance)
            rows.append(dict(case=case, families=n, replicates=replicates, target=float(truth),
                coverage=float(covered / replicates), covered=int(covered), invalid_variances=invalid,
                coverage_wilson95=wilson(covered, replicates), bias=float(np.mean(estimates) - truth),
                empirical_variance=float(np.var(estimates, ddof=1)), mean_variance=float(np.mean(variances)),
                screen_pass=bool(invalid == 0 and covered / replicates >= .925)))
    return dict(status="reference_size_sparse_screen_complete", source=identity(__file__),
        helper=identity(Path(__file__).with_name("validate_dyadic_f1_uncertainty.py")), reference=reference,
        reference_only_histogram=histogram, seed=SEED, numpy_version=np.__version__, rows=rows,
        all_cells_pass=all(r["screen_pass"] for r in rows), benchmark_intervals_admitted=False,
        limitations=["No native prediction counts or method identities are used to fit simulation parameters.",
            "Reference-size distribution only; not exact biological eligibility, annotations or evolutionary dependence.",
            "Sparse computation preserves zero contributions but cannot fix a statistically invalid Wald approximation.",
            "Rare-perfect case has fixed expected error counts as family count increases; it is a boundary stress case.",
            "All invalid variances count as noncoverage; no post-hoc clipping or correction."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = run(args.repo.resolve(), args.replicates)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")

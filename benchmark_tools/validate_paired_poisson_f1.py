"""Exact conditional Poisson boundary validation, never native F1 intervals."""

import argparse
import hashlib
import json
from numbers import Integral
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import chi2, poisson


PROTOCOL = "benchmark_tools/results/PAIRED_POISSON_F1_PROTOCOL_20261010.md"
TRUTH_MASSES = (512, 23934)
MEANS = ((0., 0., 0.), (2., 0., 0.), (1., 1., 0.), (1., 0., 1.), (.5, .3, .7), (10., 20., 12.))
ALPHA = .05
TAIL_BUDGET = 1e-12
MAX_TRIPLES = 2_000_000


def pin(path):
    path = Path(path).resolve(strict=True)
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def validate_alpha(alpha):
    if isinstance(alpha, (bool, np.bool_)) or not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("Require alpha between zero and one")


def validate_truth(truth):
    if isinstance(truth, (bool, np.bool_)) or not isinstance(truth, Integral) or not 0 < truth <= 2**52:
        raise ValueError("Require positive exact integer truth mass")


def garwood(counts, alpha):
    validate_alpha(alpha)
    counts = np.asarray(counts)
    if (not np.issubdtype(counts.dtype, np.integer) or np.any(counts < 0)
            or np.any(counts > 2**52)):
        raise ValueError("Require nonnegative exact integer counts")
    values = counts.astype(float)
    lower = np.where(values == 0, 0., chi2.ppf(alpha / 2, 2 * values) / 2)
    upper = chi2.isf(alpha / 2, 2 * (values + 1)) / 2
    if not np.isfinite(lower).all() or not np.isfinite(upper).all() or np.any(lower > upper):
        raise ValueError("Invalid Poisson quantiles")
    return lower, upper


def contrast(truth, common, left_only, right_only):
    validate_truth(truth)
    c, a, b = [np.asarray(v, dtype=float) for v in (common, left_only, right_only)]
    if any(not np.isfinite(v).all() or np.any(v < 0) for v in (c, a, b)):
        raise ValueError("Require finite nonnegative means or counts")
    mass = 2. * truth
    with np.errstate(over="ignore", invalid="ignore"):
        numerator = mass * (b - a)
        denominator = (mass + c + a) * (mass + c + b)
    if not np.isfinite(numerator).all() or not np.isfinite(denominator).all():
        raise ValueError("Contrast exceeds supported numeric range")
    return numerator / denominator


def project_box(truth, bounds):
    validate_truth(truth)
    if len(bounds) != 3:
        raise ValueError("Need three category bounds")
    for low, high in bounds:
        if np.any(np.asarray(low) > np.asarray(high)):
            raise ValueError("Reversed category bounds")
    (cl, cu), (al, au), (bl, bu) = bounds
    lower = np.minimum(contrast(truth, cl, au, bl), contrast(truth, cu, au, bl))
    upper = np.maximum(contrast(truth, cl, al, bu), contrast(truth, cu, al, bu))
    return lower, upper


def interval(truth, counts, alpha=ALPHA):
    if len(counts) != 3:
        raise ValueError("Need shared, left-only and right-only counts")
    validate_alpha(alpha)
    return project_box(truth, [garwood(value, alpha / 3) for value in counts])


def enumerate_cell(truth, means, alpha=ALPHA, tail_budget=TAIL_BUDGET):
    validate_truth(truth)
    validate_alpha(alpha)
    if len(means) != 3 or any(not np.isfinite(v) or v < 0 for v in means):
        raise ValueError("Need three finite nonnegative Poisson means")
    if not np.isfinite(tail_budget) or not 0 < tail_budget < 1:
        raise ValueError("Invalid tail budget")
    cutoffs = [int(poisson.isf(tail_budget / 3, mean)) for mean in means]
    triples = int(np.prod([cutoff + 1 for cutoff in cutoffs], dtype=object))
    if triples > MAX_TRIPLES:
        raise ValueError("Too many enumerated count triples")
    axes = [np.arange(cutoff + 1) for cutoff in cutoffs]
    bounds = [garwood(axis, alpha / 3) for axis in axes]
    pmfs = [poisson.pmf(axis, mean) for axis, mean in zip(axes, means)]
    shaped = [(lo.reshape(tuple(len(lo) if j == i else 1 for j in range(3))),
               hi.reshape(tuple(len(hi) if j == i else 1 for j in range(3))))
              for i, (lo, hi) in enumerate(bounds)]
    lower, upper = project_box(truth, shaped)
    weights = pmfs[0][:, None, None] * pmfs[1][None, :, None] * pmfs[2][None, None, :]
    target = float(contrast(truth, *means))
    covered = (lower <= target) & (target <= upper)
    widths = upper - lower
    if not np.isfinite(widths).all() or np.any(widths < 0):
        raise ValueError("Invalid projected width")
    mass = float(weights.sum())
    tails = [float(poisson.sf(cutoff, mean)) for cutoff, mean in zip(cutoffs, means)]
    omitted = float(1 - np.prod([1 - tail for tail in tails]))
    covered_mass = float(weights[covered].sum())
    rectangle_mass = float(np.prod([pmf[(lo <= mean) & (mean <= hi)].sum()
                                   for (lo, hi), pmf, mean in zip(bounds, pmfs, means)]))
    return dict(truth_mass=truth, means=list(means), alpha=alpha, target=target, cutoffs=cutoffs,
                enumerated_triples=triples, enumerated_mass=mass, omitted_tail_mass=omitted,
                covered_enumerated_mass=covered_mass, rectangle_covered_mass=rectangle_mass,
                coverage_tail_only_bounds=[covered_mass, min(1., covered_mass + omitted)],
                mean_width_enumerated_mass=float((weights * widths).sum()),
                numerical_check_pass=bool(covered_mass >= 1 - alpha - 1e-10 and omitted <= 1e-10),
                numerical_rounding_certified=False)


def shock_counterexample(truth):
    target = float(contrast(truth, 0., 10., 0.))
    outcomes = []
    for probability, left_only in ((.99, 0), (.01, 1000)):
        low, high = interval(truth, (0, left_only, 0))
        outcomes.append(dict(probability=probability, left_only_count=left_only,
                             interval=[float(low), float(high)], covers_target=bool(low <= target <= high)))
    return dict(truth_mass=truth, target=target, mean_left_only=10.,
                law="Non-Poisson common-shock mixture: zero with probability .99; 1000 with probability .01",
                outcomes=outcomes, coverage=sum(r["probability"] for r in outcomes if r["covers_target"]),
                poisson_marginal_assumption_satisfied=False, native_dependence_model=False)


def run(root, protocol_sha):
    protocol = pin(Path(root) / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Changed prospective protocol")
    rows = [enumerate_cell(truth, means) for truth in TRUTH_MASSES for means in MEANS]
    controls = [shock_counterexample(truth) for truth in TRUTH_MASSES]
    if pin(Path(root) / PROTOCOL) != protocol:
        raise ValueError("Protocol changed during enumeration")
    return dict(schema="conditional_paired_poisson_f1_v1", status="conditional_exact_boundary_enumeration_complete",
                protocol=protocol, source=pin(__file__), rows=rows, model_violation_controls=controls,
                numpy_version=np.__version__, scipy_version=scipy.__version__,
                joint_law="Independent Poisson shared/left-only/right-only whole-panel counts; fixed perfect-recall truth mass",
                target="Difference of F1 ratios at expected count vectors, not expectation of sample F1",
                all_conditional_cells_pass=all(row["numerical_check_pass"] for row in rows),
                excluded_laws=["imperfect_recall", "unequal_size_latent_node", "shared_clade", "native_VGNC", "other_QfO_endpoints"],
                overall_uncertainty_method_admitted=False, native_intervals_admitted=False, new_bootstrap_draws=0,
                new_inference_or_benchmark_scoring=False, independent_biological_confirmation=False, publication_ready=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = run(args.root, args.protocol_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(status=result["status"], cells=len(result["rows"]),
                         all_conditional_cells_pass=result["all_conditional_cells_pass"], native_intervals_admitted=False)))


if __name__ == "__main__":
    main()

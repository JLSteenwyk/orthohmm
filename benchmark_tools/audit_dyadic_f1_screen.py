"""Independent closed-form linearized variance check for the synthetic screen."""

import argparse
import hashlib
import json
from pathlib import Path


def oracle(n, case):
    if n not in (64, 256) or case not in (
            "independent_sparse", "node_sparse", "node_dense", "clade_sparse"):
        raise ValueError("Outside frozen simulation design")
    divisor = 32 if case == "node_dense" else n
    dyads = n * (n - 1) / 2
    fp = [8 * dyads / divisor, 5 * dyads / divisor]
    recall = [0.75, 0.80]
    slopes = [0., 0.] if case == "independent_sparse" else [0.15, 0.10]
    truth_pairs = 20 * n
    tp = [truth_pairs * r for r in recall]
    denominator = [truth_pairs + t + f for t, f in zip(tp, fp)]
    # Differentiate after FN = fixed truth-pair count - TP.
    a = [2 * (truth_pairs + f) / d**2 for f, d in zip(fp, denominator)]
    b = [-2 * t / d**2 for t, d in zip(tp, denominator)]
    v0 = recall[0] * (1 - recall[0]) - slopes[0]**2
    v1 = recall[1] * (1 - recall[1]) - slopes[1]**2
    covariance = recall[0] * (1 - recall[1]) - slopes[0] * slopes[1]
    diagonal_noise = truth_pairs * (a[0]**2 * v0 + a[1]**2 * v1 - 2 * a[0] * a[1] * covariance)
    poisson_noise = b[0]**2 * fp[0] + b[1]**2 * fp[1] - 2 * b[0] * b[1] * fp[1]
    effect = 0. if case == "independent_sparse" else (
        20 * (a[0] * slopes[0] - a[1] * slopes[1])
        + (n - 1) * 0.35 * (8 * b[0] - 5 * b[1]) / divisor)
    latent_variance = n * effect**2 * (8 if case == "clade_sparse" else 1)
    return dict(target=2 * tp[0] / denominator[0] - 2 * tp[1] / denominator[1],
                conditional_diagonal_variance=diagonal_noise,
                conditional_poisson_variance=poisson_noise,
                latent_effect_variance=latent_variance,
                exact_linearized_variance=diagonal_noise + poisson_noise + latent_variance)


def run(source, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    raw = source.read_bytes()
    original = json.loads(raw)
    rows = []
    for row in original["rows"]:
        expected = oracle(row["families"], row["case"])
        if abs(row["known_ratio_of_expected_counts_difference"] - expected["target"]) > 1e-14:
            raise ValueError("Synthetic target mismatch")
        rows.append(dict(case=row["case"], families=row["families"], oracle=expected,
            empirical_to_oracle_variance=row["empirical_variance"] / expected["exact_linearized_variance"],
            candidate_to_oracle_variance=row["mean_candidate_variance"] / expected["exact_linearized_variance"]))
    result = dict(source=dict(path=str(source), sha256=hashlib.sha256(raw).hexdigest()),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), rows=rows,
        status="closed_form_linearized_variance_readback", benchmark_intervals_admitted=False,
        limitations=["Post-screen algebra audit, not a new coverage gate or exact nonlinear-F1 variance.",
                     "Simulation model assumptions do not establish validity for native VGNC."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.source.resolve(), args.output.absolute())

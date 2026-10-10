"""Prospective exact-mass validation of the fixed-population native FAS design."""

import argparse
from bisect import bisect_right
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
from math import comb
from pathlib import Path

import scipy

from benchmark_tools import native_fas_sampling_interval as kernel
from benchmark_tools.probe_native_fas_context import fingerprint, require


KERNEL_SHA = "803ff0f455213e6a55b83e30c586ab1c57566830f8f877c916b78e84bf3658f9"
HELPER_SHA = "d722944a7cadeee3e25433fc472cdf724ba58309aaac29952338e5c7575bae06"
PROTOCOL_SHA = "e47b22d456dea7a41fcdf9248981d0c2d4eafec5c6dcb6ace5ad2153297b70dd"
ALPHA = Fraction(1, 20)
ERROR = ALPHA / 16
TOLERANCE = 1e-12
JOINT_INDICES = (3, 4, 5, 6, 7, 8, 9, 11)


def plan():
    half = [(2, ".25"), (2, ".75")]
    binary = [(8, "0"), (8, "1")]
    return [
        dict(name="zero_returns", pre=[(3, ".25"), (3, ".75")], missing=[(8, None)], c=4, k=3),
        dict(name="singleton_return", pre=half, missing=[(1, ".75")], c=1, k=4),
        dict(name="zero_draws", pre=half, missing=[(4, None), (4, ".5")], c=0, k=1),
        dict(name="complete_draws", pre=[(3, ".25"), (3, ".75")], missing=[(3, None), (2, ".25"), (3, ".75")], c=8, k=6),
        dict(name="random_denominator_example", pre=[(1, ".1"), (1, ".3"), (1, ".8")], missing=[(2, None), (1, ".2"), (1, ".9")], c=2, k=2),
        dict(name="balanced_partial", pre=binary, missing=[(64, None), (32, "0"), (32, "1")], c=64, k=8),
        dict(name="high_pre_low_return", pre=[(4, "0"), (16, "1")], missing=[(48, None), (64, "0"), (16, "1")], c=64, k=10),
        dict(name="low_pre_high_return", pre=[(16, "0"), (4, "1")], missing=[(48, None), (16, "0"), (64, "1")], c=64, k=10),
        dict(name="rare_partial_returns", pre=binary, missing=[(126, None), (1, "0"), (1, "1")], c=32, k=4),
        dict(name="all_return_partial", pre=binary, missing=[(64, "0"), (64, "1")], c=64, k=8),
        dict(name="native_cap_sparse", pre=half, missing=[(11998, None), (1, ".2"), (1, ".9")], c=9000, k=3),
        dict(name="native_cap_dense", pre=half, missing=[(1, None), (11998, ".25"), (1, ".75")], c=9000, k=3),
        dict(name="native_cap_no_returns", pre=half, missing=[(12000, None)], c=9000, k=3),
        dict(name="native_cap_all_returns", pre=half, missing=[(12000, ".75")], c=9000, k=3),
        dict(name="unequal_strata", pre=[(16, "0"), (48, "1")], missing=[(3, None), (1, ".25"), (1, ".75")], c=3, k=38),
    ]


def finite_outcomes(categories, draws):
    """Collapse equiprobable subsets by category counts, without a random model."""
    counts = tuple(n for n, _ in categories)
    total = sum(counts)
    require(0 <= draws <= total and all(n >= 0 for n in counts), "Invalid finite population")
    denominator = comb(total, draws)

    def visit(index, remaining, selected, multiplicity):
        if index == len(counts):
            if remaining == 0:
                returned = sum(n for n, (_, score) in zip(selected, categories) if score is not None)
                score_sum = sum((n * Fraction(score) for n, (_, score) in zip(selected, categories)
                                 if score is not None), Fraction())
                yield dict(selected=selected, returned=returned, score_sum=score_sum,
                           probability=Fraction(multiplicity, denominator))
            return
        lower = max(0, remaining - sum(counts[index + 1:]))
        for n in range(lower, min(counts[index], remaining) + 1):
            yield from visit(index + 1, remaining - n, selected + (n,), multiplicity * comb(counts[index], n))
    return list(visit(0, draws, (), 1))


@lru_cache(maxsize=None)
def exact_count_mass(M, G, c, r):
    if not max(0, c - M + G) <= r <= min(c, G):
        return Fraction()
    return Fraction(comb(G, r) * comb(M - G, c - r), comb(M, c))


@lru_cache(maxsize=None)
def exact_weight(M, G, c, k):
    return sum((exact_count_mass(M, G, c, r) * Fraction(k, k + r)
                for r in range(max(0, c - M + G), min(c, G) + 1)), Fraction())


@lru_cache(maxsize=None)
def checked_count_bounds(M, c, r):
    lo, hi = kernel.count_interval(M, c, r, float(ERROR))
    support_low, support_high = r, M - c + r
    def tail(G, upper):
        start, stop = max(0, c - M + G), min(c, G)
        return sum((exact_count_mass(M, G, c, n) for n in range(start, stop + 1)
                    if (n >= r if upper else n <= r)), Fraction())
    require(tail(lo, True) >= ERROR / 2 and tail(hi, False) >= ERROR / 2, "Rejected count endpoints")
    require(lo == support_low or tail(lo - 1, True) < ERROR / 2, "Lower endpoint not sharp")
    require(hi == support_high or tail(hi + 1, False) < ERROR / 2, "Upper endpoint not sharp")
    return lo, hi


def decimal(value):
    return Decimal(value.numerator) / Decimal(value.denominator)


def reference_interval(M, c, k, pre_mean, r, sample_mean):
    lo, hi = checked_count_bounds(M, c, r)
    if r:
        width = ((Decimal(2) / decimal(ERROR)).ln() / (2 * r)).sqrt()
        mean = decimal(sample_mean)
        mean_low, mean_high = max(Decimal(0), mean - width), min(Decimal(1), mean + width)
    else:
        mean_low, mean_high = Decimal(0), Decimal(1)
    corners = [decimal(w) * decimal(pre_mean) + (1 - decimal(w)) * mean
               for w in (exact_weight(M, lo, c, k), exact_weight(M, hi, c, k))
               for mean in (mean_low, mean_high)]
    return dict(count=[lo, hi], mean=[mean_low, mean_high], target=[min(corners), max(corners)])


def enclosed(value, bounds):
    # Reference arithmetic is 80 digits; this guards exact endpoint equalities.
    margin = Decimal("1e-60")
    point = decimal(value)
    return bounds[0] - margin <= point <= bounds[1] + margin


def validate_cell(cell):
    P = sum(n for n, _ in cell["pre"])
    M = sum(n for n, _ in cell["missing"])
    G = sum(n for n, score in cell["missing"] if score is not None)
    muP = sum((n * Fraction(score) for n, score in cell["pre"]), Fraction()) / P
    muG = sum((n * Fraction(score) for n, score in cell["missing"] if score is not None), Fraction()) / G if G else Fraction()
    c, k = cell["c"], cell["k"]
    require(1 <= k <= P and c <= min(M, 9000), "Invalid prospective requests")
    old, new = finite_outcomes(cell["pre"], k), finite_outcomes(cell["missing"], c)
    require(sum((s["probability"] for s in old), Fraction()) == 1 and
            sum((s["probability"] for s in new), Fraction()) == 1, "Finite law mass missing")
    # Independently average the actual ratio over both strata's subset outcomes.
    direct = sum((a["probability"] * b["probability"] *
                  (a["score_sum"] + b["score_sum"]) / (k + b["returned"])
                  for a in old for b in new), Fraction())
    weight = exact_weight(M, G, c, k)
    target = weight * muP + (1 - weight) * muG
    require(direct == target, "Independent native subset target disagrees")
    plugin = (k * muP + Fraction(c * G, M) * muG) / (k + Fraction(c * G, M))
    count_covered = mean_covered = target_covered = Fraction()
    max_error = abs(kernel.expected_native_mean(M, G, c, k, float(muP), float(muG)) - float(target))
    states = []
    for state in new:
        r = state["returned"]
        mean = state["score_sum"] / r if r else None
        actual = kernel.method_interval(M, c, k, float(muP), r, float(mean) if r else None, float(ERROR))
        oracle = reference_interval(M, c, k, muP, r, mean)
        require(actual["success_count_bounds"] == oracle["count"], "Count projection differs")
        for key, reference in (("return_mean_bounds", oracle["mean"]), ("expected_native_mean_bounds", oracle["target"])):
            max_error = max(max_error, *(abs(x - float(y)) for x, y in zip(actual[key], reference)))
        prob = state["probability"]
        count_covered += prob if oracle["count"][0] <= G <= oracle["count"][1] else 0
        mean_covered += prob if not G or enclosed(muG, oracle["mean"]) else 0
        covered = enclosed(target, oracle["target"])
        target_covered += prob if covered else 0
        states.append(dict(probability=prob, interval=actual, reference=oracle["target"], covered=covered))
    require(count_covered >= 1 - ERROR and mean_covered >= 1 - ERROR,
            "Component confidence coverage failed")
    require(target_covered >= 1 - 2 * ERROR and max_error <= TOLERANCE, "Target coverage/numerics failed")
    return dict(name=cell["name"], population=dict(P=P,M=M,G=G,c=c,k=k),
                literal_native_cap=c == min(M,9000), subset_categories=len(new),
                two_stratum_outcomes=len(old)*len(new), target_fraction=str(target),
                native_ratio_average_fraction=str(direct), expected_count_plugin_fraction=str(plugin),
                plugin_minus_target_fraction=str(plugin-target), count_coverage_fraction=str(count_covered),
                mean_coverage_fraction=str(mean_covered), target_coverage_fraction=str(target_covered),
                numerical_max_error=max_error), states, target


def validate_joint(rows, distributions, targets):
    selected = [distributions[i][::(-1 if j % 2 else 1)] for j, i in enumerate(JOINT_INDICES)]
    cuts, accumulated = {Fraction(0),Fraction(1)}, []
    for states in selected:
        total, endpoints = Fraction(), []
        for state in states:
            total += state["probability"]
            endpoints.append(total)
        require(total == 1, "Invalid coupled marginal")
        accumulated.append(endpoints)
        cuts.update(endpoints)
    cuts = sorted(cuts)
    means_covered = differences_covered = Fraction()
    max_error = 0.
    for lower, upper in zip(cuts,cuts[1:]):
        middle = (lower+upper)/2
        states = [dist[bisect_right(cdf,middle)] for dist,cdf in zip(selected,accumulated)]
        means_ok = all(s["covered"] for s in states)
        pairs_ok = True
        for a,b in combinations(range(8),2):
            oracle = [states[a]["reference"][0]-states[b]["reference"][1],
                      states[a]["reference"][1]-states[b]["reference"][0]]
            actual = kernel.difference_interval(states[a]["interval"],states[b]["interval"])
            reverse = kernel.difference_interval(states[b]["interval"],states[a]["interval"])
            require(reverse == [-actual[1],-actual[0]], "Method swapping failed")
            max_error=max(max_error,*(abs(x-float(y)) for x,y in zip(actual,oracle)))
            pairs_ok = pairs_ok and enclosed(targets[JOINT_INDICES[a]]-targets[JOINT_INDICES[b]],oracle)
        means_covered += upper-lower if means_ok else 0
        differences_covered += upper-lower if pairs_ok else 0
    union_failure = sum((2-Fraction(rows[i]["count_coverage_fraction"])-
                         Fraction(rows[i]["mean_coverage_fraction"]) for i in JOINT_INDICES),Fraction())
    require(union_failure <= ALPHA and means_covered >= 1-ALPHA and
            differences_covered >= 1-ALPHA and max_error <= TOLERANCE, "Joint coverage/projection failed")
    return dict(method_count=8,differences=28,rank_coupling_segments=len(cuts)-1,
                coupled_mean_coverage_fraction=str(means_covered),coupled_difference_coverage_fraction=str(differences_covered),
                arbitrary_dependence_union_failure_bound_fraction=str(union_failure),numerical_max_error=max_error,
                scope="Common-uniform rank coupling with alternate reversed outcome order; uniform precomputed selections independent within each method.")


def context_control():
    # A returns only when B accompanies it; no fixed per-pair return set exists.
    scores={"A":Fraction(1,5),"B":Fraction(4,5),"C":None}
    values=[]
    for subset in combinations(scores,2):
        returned=[scores[x] for x in subset if scores[x] is not None and (x!="A" or "B" in subset)]
        values.append((Fraction(2,5)+sum(returned))/ (1+len(returned)))
    direct=sum(values)/len(values)
    singleton_target=exact_weight(3,1,2,1)*Fraction(2,5)+(1-exact_weight(3,1,2,1))*Fraction(4,5)
    require(direct==Fraction(22,45) and singleton_target==Fraction(8,15) and direct!=singleton_target,"Control failed to violate fixed returns")
    return dict(status="fixed_return_assumption_rejected", actual_expected_ratio_fraction=str(direct),
                singleton_fixed_return_target_fraction=str(singleton_target),
                coverage_claimed=False,native_admission=False)


def run(protocol):
    source=fingerprint(kernel.__file__)
    binding=fingerprint(protocol)
    helper=fingerprint(fingerprint.__code__.co_filename)
    require(source["sha256"]==KERNEL_SHA and binding["sha256"]==PROTOCOL_SHA and
            helper["sha256"]==HELPER_SHA,"Changed frozen kernel/protocol/helper")
    rows,distributions,targets=[],[],[]
    with localcontext() as ctx:
        ctx.prec=80
        for cell in plan():
            row,states,target=validate_cell(cell)
            rows.append(row);distributions.append(states);targets.append(target)
            print(cell["name"],row["target_coverage_fraction"],flush=True)
        joint=validate_joint(rows,distributions,targets)
        control=context_control()
    require(fingerprint(kernel.__file__)==source and fingerprint(protocol)==binding and
            fingerprint(helper["path"])==helper,"Inputs changed during validation")
    return dict(schema="native_fas_finite_sampling_validation_v1",status="prospective_validation_complete",
                kernel=source,protocol=binding,provenance_helper=helper,driver=fingerprint(__file__),plan=plan(),alpha_fraction=str(ALPHA),
                component_error_fraction=str(ERROR),decimal_precision=80,reference_endpoint_margin="1e-60",
                numerical_tolerance=TOLERANCE,scipy_version=scipy.__version__,cells=rows,joint=joint,context_control=control,
                sampling_probabilities="Exact rational subset multiplicities; no simulation, omitted tail or new RNG draw.",
                coverage_scope="Exact rational covered-outcome mass evaluated with 80-digit reference endpoints; floating kernels compared separately.",
                historical_scores_rerun=False,native_sampling_law_admitted=False,native_intervals_admitted=False,publication_ready=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result=run(args.protocol)
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(status=result["status"],cells=len(result["cells"]),joint=result["joint"])))


if __name__=="__main__":
    main()

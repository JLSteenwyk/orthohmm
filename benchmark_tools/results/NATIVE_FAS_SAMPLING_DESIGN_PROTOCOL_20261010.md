# Native FAS Sampling Design Investigation

Prospective work for original requirement 4.1. No estimator has run and no native
interval is admitted. The aim is uncertainty for the expected native
post-attrition FAS statistic under its actual sampling algorithm, not a new
complete-population or shared-pair endpoint.

## Inspected Native Code

The retained benchmark scorer shuffles its precomputed pair list, requests
k pairs using its rounded native calculation, and separately shuffles the
lookup-missing list to request c=min(M,9000) pairs. The final mean includes k
precomputed scores and the r newly returned numeric scores. Missing numeric
returns are silently omitted. The function's locally sliced missing list
does not replace the caller's full list; the final loop visits that full list
but includes only lookup-present pairs.

Static reading inside the retained container confirms the historical source
identities recorded by the earlier omission probe:

| Source | SHA256 |
| --- | --- |
| fas_benchmark.py | 1045c57f4d0f4787bec3d1f0690799df63c68dcc67ccfd5a925c338e3d33661d |
| calcFASmulti.py | 7f9621c0dc6737c8b928c0cf35575101decc4af1cf039e4f3c044e80ce460b75 |
| calcFAS.py | 40a920fcce24d76e7441f5f9368aa6232d779815956742c725ca16a169980c2b |
| mainFAS/greedyFAS.py | f82a78d2787d30f106a58c760cae49d28f28d55956ad188b0b259ec8b9f0f7e8 |

The child groups selected pairs by taxa, but loads complete annotation files.
Native options set MS_uni=1 without a reference proteome, so the branch that
updates reference feature counts is inactive. Pair preparation and path-limit
rejection use a protein's annotation and fixed options. Priority mode is
enabled by default and the command does not disable it; the inspected
timeout-to-priority branch is inside the non-priority branch. This is source
evidence, not a claim that every downstream helper or historical failure has
been proved sample-context invariant.

The [earlier probe](QFO_FAS_OMISSION_MECHANISM_20260928.md) establishes one
per-protein path-limit omission mechanism, not attribution of historical
omissions. That probe must not be repeated. The
[completed population audit](QFO_FAS_POPULATION_COMPLETED_22383.md) and
[stratum arithmetic](QFO_FAS_STRATUM_WEIGHT_RESULT_20261001.md) already supply
retained counts and sums; do not repeat their expensive joins.

## Required Design Assumptions

Condition on the fixed predicted pair sets, annotation files, lookup and scoring
options. Require uniform sampling without replacement within the two strata
and independent stratum selections within a method. Historical RNG states are
unavailable: this is a probabilistic design assumption supported by the shuffle
implementation, not exact replay or a proof of ideal PRNG randomness.

Additionally require that, on normal successful execution, each lookup-missing
pair has a fixed numeric-return indicator and, when returned, a fixed score
in [0,1] independent of its sampled companions, order or worker assignment.
Missingness may depend arbitrarily on protein features and scores. It must not
depend on which other pairs were sampled. Batch crashes, partial file writes,
shared mutable scoring state or context-sensitive score changes violate this
assumption and are not absorbed into a hypergeometric return model.

Do not assume biological independence of gene pairs or families. Their scores
are fixed in this design. No conclusion about sampling new genes, families or
clades follows. Cross-method randomization need not be independent for a
simultaneous confidence-set construction using a union bound.

## Exact Native Statistic Target

Let P and M be the native eligible precomputed and lookup-missing counts.
Let k be the actual rounded request count, with 1<=k<=P, and c=min(M,9000).
All retained methods have positive k. Let mu_P be the full precomputed mean
from the retained population sum, not its realized sampled mean.

Among M missing-lookup pairs, let G have fixed numeric returns and let mu_G
be their mean score. Under the stated design, R has the hypergeometric law
for c draws from M items containing G successes. Conditional on R=r>0, the
returned pairs form a uniform r-subset of those G pairs.

**Project derivation:** for the actual reported statistic Z,

```text
Z = (sum(precomputed sample) + sum(returned new sample)) / (k + R)
E[Z | R] = [k*mu_P + R*mu_G] / (k + R)
theta(G,mu_G) = w(G)*mu_P + [1-w(G)]*mu_G
w(G) = sum_r Pr(R=r | M,G,c) * k/(k+r).
```

The expectation retains the random denominator and mixture. It is NOT obtained
by plugging E[R] into that denominator, weighting by P/(P+M), completing NA
values, conditioning on shared sampled pairs or relabeling a diagnostic as a
native endpoint. G=0 gives w=1 and theta=mu_P; mu_G is then irrelevant.

The target is the expected repeated native sampled statistic under the fixed
design and successful fixed-pair scorer, not the particular fixed observed
mean or general biological accuracy. Observed published means remain unchanged.

One rational-arithmetic derivation check enumerates all 18 subsets with
precomputed scores (0.1,0.3,0.8), two requested precomputed pairs, missing-lookup
outcomes (NA,0.2,0.9,NA) and two requested new pairs. Direct averaging of the
native ratio gives 107/240, matching the hypergeometric-mixture expression;
substituting the expected return count gives 9/20 instead. This is a checked
algebraic example, not a coverage panel or native model validation.

## Candidate Confidence Construction

Invert two-sided hypergeometric count tails to form a confidence set for G.
[SciPy 1.15.3 hypergeom](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.stats.hypergeom.html)
uses total population M, success count G and draw count c. Check integer
support r<=G<=M-c+r and both exact boundary tails before any native application.

For r>0, bound mu_G from the returned score mean by sampling-without-replacement
concentration. [Bardenet and Maillard, Proposition 1.2](https://arxiv.org/pdf/1309.4029v2)
gives a Hoeffding bound for a uniform sample from fixed bounded values. Using
both tails yields half-width sqrt(log(2/delta)/(2*r)) for [0,1] scores.
For r=0 retain the entire [0,1] range, not a zero-width plug-in bound.
This conditional-on-r guarantee also bounds unconditional error by averaging
over r; use no missing-at-random assumption.
The single-valid-pair population is a degenerate case: a positive return
observes its mean exactly. Do not invoke a theorem requiring population size
at least two for that boundary.

Allocate total alpha=0.05 across the two uncertain components for each of the
eight retained methods: delta=alpha/16 per count or score-mean interval.
A union bound covers the full parameter rectangle even with arbitrary
cross-method dependence. k and the full precomputed mean are fixed here.

**Project projection:** increasing G stochastically increases R and decreases
w(G). Theta increases with mu_G; its direction in G depends on mu_G-mu_P.
Evaluate the four endpoint combinations to project each rectangle. Joint
method-mean coverage then supplies all 28 difference intervals simultaneously
by lower_left-upper_right and upper_left-lower_right. This derivation must be
tested independently; it is not an executed native result.

## Validation Before Use

First implement and test a new order/companion/worker invariance check using
the actual native scoring functions and small controlled annotations, including
numeric and path-limit-NA cases. It must differ from, and not rerun, the
completed omission probe. Inspect remaining score helpers for shared state,
sample-dependent inputs or reachable timeout changes. Do not infer invariance
from a shuffle call or one per-protein rejection branch alone.

If that prerequisite supports this direction, commit and push a tested
prospective implementation before any new numerical/native application.
Validate count-tail inversion and projected targets against exact finite
enumeration of small populations, including zero returns, complete returns,
score-dependent omissions, unequal stratum sizes and method-swapping.
Independently sum the native ratio statistic over actual subsets; reject
expected-count substitution. Retain an explicitly context-dependent return
control outside the design, rather than calling it a covered regime.

Freeze any native reporting protocol before viewing intervals. Reuse the
retained eight-method population/sample reports and their existing bindings.
No FAS rescore, new annotations, new RNG sample, historical retry or expensive
population recount is required. Preserve every observed native mean and original
endpoint. If source invariance, data mapping or target assumptions cannot be
justified, retain that specific limitation instead of executing another
convenient restricted model. No TreeFam/VGNC/GO/EC or secondary-mean uncertainty,
independent biological confirmation or publication readiness follows from this
FAS-only route.

# Exact Paired F1 Boundary Construction

Prospective mathematical validation for original goal 4.1. This does not replace
the failed dyadic screens or supply native benchmark intervals. Do not fit this
model to native errors or call its conditional target biological generalization.

## Motivation And Source Check

The previously unavailable full network paper is now accessible through
[LSE STICERD](https://sticerd.lse.ac.uk/dps/em/em629.pdf). Assumption A, Theorem 2
and Sections 2.2--2.3 require more than sparse storage; the severely sparse
regime has positive limiting average degree. The vector/smooth-function results
also need nonsingular covariance and regular derivatives. Their applicability
to the project's weighted, paired count statistic is not established.

**Project deduction:** a fixed expected number of distinct error links has
edge probability O(n^-2), hence average degree tending to zero. That boundary
is outside the paper's positive-limit regime. It is not repaired just by naming
modified jackknife empirical likelihood. The previous rare screen and the
shared-clade failure remain retained evidence, not experiments to repeat here.

## Conditional Joint Law And Target

Fix positive integer truth mass T. Both hypothetical methods recover all T
true pairs and have zero false negatives. Observe shared false-positive count
C and method-only counts A and B. Assume independent Poisson category counts
with respective means c, a, b. Counts are whole-panel random quantities: this
is NOT an assertion that native gene pairs or reference families are IID.

FP_left=C+A; FP_right=C+B; TP_left=TP_right=T; FN_left=FN_right=0.
The target is the difference of F1 ratios evaluated at expected count vectors:

```text
delta(T,c,a,b) = 2T/(2T+c+a) - 2T/(2T+c+b)
              = 2T(b-a)/[(2T+c+a)(2T+c+b)].
```

It is not E[F1_left-F1_right], uncertainty about the observed fixed-table score,
an interval conditional on selected errors, or a general native F1 model.

## Construction And Coverage Argument

Use equal-tailed [Garwood Poisson limits](https://jhanley.biostat.mcgill.ca/bios601/Intensity-Rate/Garwood1936.pdf)
for each mean, with total error allocation alpha/3 per category. The
[PDG Statistics review, equations 40.77a--b](https://pdg.lbl.gov/2021/reviews/rpp2021-rev-statistics.pdf)
provides the chi-square quantile form. Zero observations have lower mean zero
and a positive upper bound, never a plug-in zero-width interval.

**Project derivation:** the union bound makes the three-mean rectangle cover
with probability at least 1-alpha. This simultaneous-coverage step needs valid
Poisson marginals, not category independence; independence is assumed for the
joint enumeration below. Project that rectangle through the actual ratio
difference, not a linearized score or a difference of separately sampled methods.

Delta decreases with a and increases with b. Its c derivative is
2T[(2T+c+b)^-2-(2T+c+a)^-2], so extrema in c are at its endpoints (or constant
when a=b). Thus the sharp rectangle lower bound uses a_upper/b_lower and the
smaller of its two c endpoints; the upper bound uses a_lower/b_upper and the
larger endpoint value. Inclusion of the true rectangle guarantees inclusion of
the true delta. Conservatism is accepted; no post-result narrowing is allowed.

## Frozen Numerical Panel

Enumerate the independent Poisson probability mass, without random draws:

- Truth masses: 512 and 23,934. The latter is the already-retained asserted
  reference-pair mass, not an error-rate fit; it does not make the model native.
- Means (c,a,b): (0,0,0), (2,0,0), (1,1,0), (1,0,1), (0.5,0.3,0.7), (10,20,12).
- Alpha=0.05; omit at most 1e-12 probability mass using per-axis Poisson tails.
- Maximum 2,000,000 enumerated count triples per cell; one numeric thread.
- Retain all 12 cells, target, count cutoffs, enumerated mass, omitted tail,
  covered mass, widths, and rectangle-coverage mass. Tail-only coverage bounds
  do not certify floating-point rounding. No Monte Carlo confidence interval.
- A necessary numerical check is covered enumerated mass >= 0.95-1e-10 and
  omitted tail <= 1e-10, with finite nonnegative widths. Failed cells stay failed.
- Component tail identities, sharp-box endpoints, zero observations, method
  swapping and direct scalar mass summation are tested separately.

Retain two additional explicit model-violation controls, one per truth mass:
C=B=0; A=0 with probability 0.99 and A=1000 with probability 0.01. The same
mean-count target has a=10, but this marginal is not Poisson. Evaluate both
outcomes and sum their exact coverage probabilities. This panel-level common
shock represents a dependence/overdispersion warning, not a fitted biological
clade model. Report the control even if coverage is zero; never relabel it a
covered regime or remove it from the recorded validation.

Regular imperfect-recall, unequal-size latent-node and shared-clade laws from
earlier validation remain outside this model. They are explicitly NOT covered
or counted as successes. No overall uncertainty-method admission follows, even
if all conditional cells pass. Native GO/EC/FAS/TreeFam and the secondary mean
are not transformed into this model or replaced with conditional endpoints.

Commit/push tested source and this protocol before one fresh enumeration output.
Bind the explicit protocol digest and source bytes; refuse occupied outputs.
No inference, timing, score, bootstrap draw, parameter tuning, TreeFam search,
DGX access or private platform change is authorized by this numerical panel.

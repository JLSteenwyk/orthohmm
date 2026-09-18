# VGNC Dependency Structure

## Executed Audit

The reference contains16863 family labels,36986 proteins and23934 asserted
pairs. Eleven proteins have multiple family labels. Merging only labels that
share a reference protein yields16844 blocks, with11 multi-label components.
This construction uses the reference, not prediction outcomes. All scored
TP/FN rows fall within these blocks in the four audited historical stages.

| Historical stage | Within-block FP | Cross-block FP | Distinct cross-block links | Blocks incident to links | Largest prediction-link component |
|---|---:|---:|---:|---:|---:|
| multipass | 2 | 120804 | 56418 | 9445 | 371 |
| multipass_refined | 2 | 15788 | 7614 | 7215 | 13 |
| strict_profiles | 2 | 121468 | 56855 | 9470 | 374 |
| strict_profiles_refined | 2 | 15804 | 7615 | 7235 | 11 |

Counts exactly reconstruct the retained TP/FP/FN inventory. Inputs were
checked by SHA-256 before/after the audit. Machine-readable evidence is
vgnc_family_dependencies_20260918.json; executable implementation is
audit_vgnc_family_dependencies.py. Eight focused tests cover transitive
reference overlap, isolated blocks, cross-block links and malformed evidence.
This supplements the earlier exact prediction rescore, not a new database
rescore or corrected-release accuracy evaluation.

## Statistical Consequences

Merging the eleven overlapping protein memberships does not make false
positives belong to independent single-family units. Assigning each cross-block
FP to one endpoint arbitrarily, or dropping it, would not reproduce the native
precision denominator. Counting it twice would also change the statistic.
The block-pair dependency must be addressed explicitly in a paired analysis.

The prediction-link connected components above are diagnostics, NOT proposed
independent units: their membership depends on the method and outcomes, and
their sizes change substantially across stages. Reference blocks themselves
also need an exchangeability justification; disjoint proteins do not prove
independent evolution or prediction errors.

[Owen's pigeonhole bootstrap](https://arxiv.org/abs/0712.1111) treats crossed
row/column dependence under specified models and cautions against naive
resampling. It does not by itself validate our symmetric same-population
family-pair setting, diagonal contributions, nonlinear F1, or outcome-defined
blocks. Applying it unchanged would be an unsupported inference here.

Next statistical work must define the sampling target, handle within-block
and cross-block contributions consistently, preserve shared resampling across
methods, recompute the native ratio statistics and validate coverage under
relevant dependence models before claiming confidence intervals. Until then,
VGNC remains a point estimate with a documented uncertainty limitation; this
does not license omission of negative results or imply that paired inference
is impossible. No scores or scientific method settings were changed.

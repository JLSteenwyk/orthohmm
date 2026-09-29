# All-Method OrthoBench Descriptive Extension

Before joining the five additional comparator methods to the frozen feature
bins, retain this protocol in a pushed commit. Use the existing fourteen
family strata and definitions unchanged: family size, copy number, alignment
identity, relative length and composition. No bin boundary or reference
membership may be selected using the new method outcomes.

Include all eight retained methods: OrthoHMM high sensitivity and satellite_v2,
full and sequence-only OrthoFinder 3.1.5, SonicParanoid 2.0.9, ProteinOrtho
6.3.6, FastOMA 0.3.5 and OrthoMCL 1.4. Source the original three methods from
the pinned paired-uncertainty report and the other five from the pinned
retained-comparator readback. Require identical RefOG names/sizes across them.
Preserve full-reference per-family TP/FP/FN and low-certainty conventions.

Within each method/stratum, sum the existing family-size-weighted sufficient
statistics and recompute precision, recall and F1. This is not an arithmetic
mean of family F1. Export all 112 method/bin rows, including empty strata as
missing rather than zero. Report family counts and preserve weighted counts.
Require reproduction of every original three-method point estimate to 1e-10
percentage points. Recheck source and prediction-file identities; no native
inference or scorer rerun is needed for this descriptive export.

Do not calculate new intervals, significance tests, ranks or selected winning
strata. The original three-method paired intervals and their 84-endpoint
multiplicity treatment remain unchanged and do not apply to the five added
comparators. These are development-exposed descriptive results, not independent
confirmation, mechanism evidence or a basis for tuning defaults. Verify the
export arithmetically and retain all method rows before summarizing it.

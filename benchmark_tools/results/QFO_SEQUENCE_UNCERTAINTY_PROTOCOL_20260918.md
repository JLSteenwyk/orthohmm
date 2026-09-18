# Corrected QfO Sequence-Control Uncertainty

## Scope And Freeze

This protocol extends SEQUENCE_SEARCH_CONTROL_PROTOCOL_20260916.md and
QFO_SEQUENCE_SEARCH_CONTROL_20260918.md before corrected sequence-control
graphs or scores exist. OrthoBench sequence-control results and historical
QfO outcomes have already been inspected. This is development-exposed,
exploratory evidence, not independent validation or selection-adjusted inference.
It does not authorize parameter tuning or a new publication-method default.

The three predictions, in fixed order, are corrected-release initial HMM
`p0_c0_r0`, DIAMOND `all_hits`, and DIAMOND `top100`. Profile expansion,
candidate expansion and phylogenetic reconciliation are off. Both DIAMOND
variants retain the frozen downstream graph/refinement settings. All-hit
search is the principal control; post-search top100 is diagnostic, not a
simulation of the HMM prefilter. Equal E-values do not establish equal
sensitivity, calibration or computational effort.

## Admission And Statistic

Require completed, independently admitted corrected-release assessments for
all three predictions. Bind admission, execution, converted predictions and
native SwissTrees raw/summary files by path and SHA-256. Reject incomplete
panels, historical prediction counts, mismatched variants and changed inputs.
The count assembler must independently reconstruct each admitted SwissTrees
endpoint within absolute tolerance 1e-6 and recheck source hashes afterward.
Historical evidence may anchor reference identities, never prediction counts.

Use the frozen SwissTrees reference SHA-256
`40528e3537ba57e345c443c8024f99dc3940a66e2d102703754b41bb29858ecd`.
Require all 18 reference families, the same represented genes and truth labels
across methods, 10,765 reference relations and no overlapping represented
genes. Any contradiction stops this analysis; do not silently drop families.

For each family, divide each raw confusion count by two and add one, matching
the native scorer. Calculate precision TP/(TP+FP) and recall TP/(TP+FN).
Average these quantities equally across families; aggregate F1 is the harmonic
mean of macro precision and macro recall, NOT mean family F1 or pooled-pair F1.

## Paired Resampling

Draw 18 families with replacement in each of 100,000 replicates using NumPy
PCG64, seed 20260923, with multinomial family multiplicities. Share every draw
across all three predictions and recompute the actual aggregate statistic.
Record NumPy version, seed, implementation and protocol hashes.

Report exactly two oriented contrasts: all_hits minus p0_c0_r0, and top100
minus p0_c0_r0. For each report F1, precision (PPV), and recall (TPR) differences
in raw 0-to-1 units, nominal 95% percentile intervals and Bonferroni intervals
over all six endpoints. Adjusted quantiles are 0.05/12 and 1-0.05/12; use linear
quantile interpolation. Do not select the best variant or add a top100 versus
all_hits contrast after inspecting results.

Also report per-family differences and wins/ties/losses (absolute tolerance
1e-10). These are descriptive; family F1 differences do not sum to the
aggregate F1 difference. Zero-prediction outcomes retain the native prior;
failed scoring is missing evidence, not an imputed zero or a dropped family.

## Interpretation And Outstanding Work

Intervals are conditional on 18 curated reference families. Disjoint genes
do not prove family exchangeability: shared evolution and merged predictions
can induce dependence. Adjustment covers these six endpoints only, not the
full history of development, other QfO challenges or the custom six-metric
mean. Do not manufacture TreeFam family units from pooled components.

Report all six native QfO endpoints and coverage descriptively when admitted.
GO, EC, VGNC, TreeFam and FAS need separately defensible uncertainty methods;
SwissTrees intervals do not satisfy that remaining requirement. Search overlap,
hit counts, threshold distributions and resource evidence are separate results.
Concurrent shared-host runs cannot establish matched-efficiency superiority.

The numerical engine is preparatory only. A source-bound count assembler and
production runner remain required before any computed interval is admitted.

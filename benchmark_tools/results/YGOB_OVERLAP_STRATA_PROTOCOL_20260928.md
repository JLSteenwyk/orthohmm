# YGOB Overlap-Stratified Diagnostic

Freeze this specification before computing the new subgroup scores. Overall
YGOB results and the historical overlap counts have already been inspected.
This is a secondary descriptive error analysis, not a new independent test,
a replacement endpoint, a prespecified primary analysis, or a basis for tuning.

## Inputs And Fixed Groups

- Retained score summary: `ygob_frozen_results_20260916.json`, SHA-256
  `3927d81b1fd851eef3c8ce1343679558ef4f93bcfe2e643435e369587a44a4ef`.
- Full sufficient statistics: the exact `full_results` reference in that
  summary, SHA-256
  `5104e000d0c3bb8687103646be8011fbc97f8b05f7eb0a925bbfdd1e76200e83`.
- Screen: `ygob_homology_screen_20260916.json`, SHA-256
  `c142f6ecba369c8192eea131de0827eafabd498399ad25c80b7136b33d34a74f`.
- Overlap admission: `ygob_overlap_admission_20260916.json`, SHA-256
  `d2a8aa5d2901c15e50a162aae8c53c54df261753d5a4c3999e89e641477f7054`.

Use the existing screen-positive pillar IDs exactly as retained. Screen-negative
means their complement within the scored reference universe. Validate unique
IDs and membership; reject missing or unexpected pillars. Do not change search
thresholds, screen labels, protein inputs or inference predictions. The existing
screen used DIAMOND very-sensitive, E <= 1e-5, identity >= 30%, and both query
and target coverage >= 50%, retaining one best qualifying hit per query.
A negative result does not establish lack of remote homology or independent
annotation ancestry. No new accuracy-based filtering is permitted.

Report all four retained methods: OrthoHMM high sensitivity, OrthoHMM
satellite_v2, full OrthoFinder 3.1.5 and its sequence-only checkpoint. Preserve
the original root-level group-recovery interpretation, including within-species
co-membership, rather than calling it resolved pairwise orthology.

## Statistics And Checks

Partition the original per-pillar TP/FP/FN sufficient statistics, not the
prediction files or gene universe. Retain each pillar's original allocated FP
count, including half of every cross-pillar pair touching another stratum.
Filtering predictions before scoring would remove those penalties and change
the endpoint. Recompute F1, precision and recall from each stratum's summed
counts; do not average pillar F1. Report sizes, reference genes, truth-pair
counts, gene coverage and exact-pillar counts alongside the scores. Include
zero-truth-pair pillar counts so undefined per-pillar ratios remain visible.

Require that strata counts, genes, coverage and exact-pillar counts sum to the
full retained totals for every method, and that reconstructed overall metrics
match the admitted summary. Validate finite nonnegative integer/half-integer
counts, TP+FN = choose(genes,2), valid coverage, unique pillars and identical
reference signatures across methods. Preserve the frozen zero-denominator
convention (score zero), but label this convention rather than treating zero
as measured accuracy when the ratio is undefined.

Report descriptive percentage-point differences against full OrthoFinder for
both strata. Do not compute new confidence intervals or significance tests:
cross-pillar errors and common ancestry do not become independent because
labels were partitioned. Do not interpret differences causally as effects of
development exposure, evolutionary distance or homology sensitivity.

## Interpretation And Deliverables

Produce tested aggregation code, a machine-readable provenance-bound result,
and a table covering both strata and all methods, retaining negative findings.
This can identify where the existing recall/precision trade-off occurs. It
cannot convert development-exposed results into independent-family validation.
No method setting, primary score, historical artifact or superiority claim
will be changed based on this diagnostic. Any future method change still
requires separate confirmation under the publication goal.

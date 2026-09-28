# YGOB Overlap-Stratified Diagnostic

The [specification](YGOB_OVERLAP_STRATA_PROTOCOL_20260928.md) was committed and
pushed as e51ffc54 before these subgroup scores were computed. Overall YGOB
results and overlap counts were already known, so this remains a secondary
descriptive analysis, not a fresh independent confirmation. No inference was
rerun and no method settings, reference universe or screen thresholds changed.

The [generated all-method table](YGOB_OVERLAP_STRATA_TABLE_20260928.md) and
[machine-readable result](ygob_overlap_strata_20260928.json) retain both strata.
The existing overlap verifier freshly checked the original development inputs,
query inputs, reference, hit table, source/executable identities and scheduler
record; its admission reproduced exactly. This is retained-hit verification,
not a new search. Full per-pillar statistics match the pinned admitted summary;
both strata reconstruct every method's overall counts and coverage exactly.

## Findings

The table reports ratios of summed retained counts, not mean pillar F1 or scores
from rescoring a smaller gene universe. The generated detailed table contains
precision, recall, coverage and exact-pillar counts for every row.

Phylogenetic OrthoHMM minus full OrthoFinder is +12.916292 F1 percentage points
in the screen-negative stratum, but its recall is 23.929262 points lower
(60.099475% versus 84.028737%). Its precision is 55.916909% versus 30.743257%.
In the screen-positive stratum, its F1 is 1.634631 points lower and recall is
6.382266 points lower than full OrthoFinder. Thus the higher screen-negative
F1 is accompanied by a precision/recall trade-off; it is not evidence of
superior recovery of remote orthologs.

The strata are compositionally very different. Screen-negative has 3,298
pillars but only 6,104 reference genes and 18,095 truth pairs; 2,893 pillars
are singletons. Screen-positive has 6,952 pillars, 77,287 reference genes and
578,995 truth pairs; 2,017 pillars are singletons. Singleton truth contributes
no true pair, but erroneous merging can still contribute allocated false
positives. The difference in family-size composition prevents a simple causal
interpretation in terms of development exposure or evolutionary distance.

All four methods retain 100% reference-gene coverage in both strata. Presence
in some group does not imply correct co-membership. Original cross-pillar
false positives remain allocated half to each endpoint pillar, including
cross-stratum pairs. No predictions were filtered to remove those penalties.

## Validation And Interpretation

33 focused aggregator/scorer tests passed, including cross-stratum penalties,
half-counts, zero-denominator flags, empty strata, incorrect totals/signatures,
duplicate labels, invalid counts and order invariance. A separate
[exact-rational arithmetic check](ygob_overlap_strata_arithmetic_20260928.json)
recomputed all 24 reported metric values from TP/FP/FN and matched within
1e-15. It checks arithmetic, not native output reconstruction or independence.

No confidence interval or significance claim is made. Screen-negative does
not mean family-disjoint: remote homology, shared curation and other exposure
remain possible. This diagnostic supplies additional error characterization
for the frozen method; it does not repair the independent-family evidence
gap, revise the overall YGOB ranking or authorize outcome-driven tuning.
Neither this result nor the completed collector fixtures establishes
publication readiness.

Reproduce aggregation in a fresh output file:

```bash
python -B -m benchmark_tools.summarize_ygob_overlap_strata \
  --root . --output /fresh/ygob-overlap-strata.json
```

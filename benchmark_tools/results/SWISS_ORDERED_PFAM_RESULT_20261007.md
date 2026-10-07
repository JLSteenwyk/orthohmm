# Ordered-Pfam Annotation Error Analysis

The prospective [protocol](SWISS_ORDERED_PFAM_PROTOCOL_20261007.md) was pushed
as 6948fb7b before feature construction. Both sources and 56 invented tests
were pushed as d09e2bc41e12020c6217e28e22811c67477eaee4 before the selected
invocations. Joined descriptor/existing-rational tests: 74 passed in 1.38 s.
No method change, annotation extraction, alignment, tree inference, raw scorer,
bootstrap or new OrthoFinder run occurred.

## Retained Annotation Descriptors

All 18 corrected SwissTrees families and all 563 canonical proteins are retained.
Every annotation length equals its retained ungapped protein length. A separate
all-pair interval check agrees with the primary adjacent-interval rule for every
protein: 434 have usable ordered signatures, 124 have ambiguous coordinates and
five have zero Pfam hits. Ambiguous and zero-hit states are not biological loss
or experimentally verified incomplete architecture.

The fixed partition is four same-signature families (CITE, POP, RPS, SUMF),
one all-usable multiple-signature family (Clusterin), and 13 families with at
least one unusable member. Among all 5,564 within-family usable-protein pairs
with identical domain type/multiplicity multisets, none differs in order.
These dependent pairs are not independent observations for confidence intervals.
Clusterin's two signatures are one versus two `pfam_Clusterin` instances, not a
rearrangement of a shared multiset. Other different signatures can likewise
reflect content/repeat differences, missing hits and annotation assumptions.
No negative finding proves biological conservation or absence of rearrangement.

## Descriptive Conditional Effects

The [complete machine-generated table](native_swiss_ordered_pfam_20261007_v1/TABLE.md)
retains all 12 scores and eight contrasts. The original TP/2+1, FP/2+1, FN/2+1
prior is applied per family; precision and recall are averaged across families
before harmonic F1. Counts are not pooled and family F1 is not averaged.

| Annotation-defined bin | Families | R effect on F1 (pp) | C effect on F1 (pp) |
| --- | ---: | ---: | ---: |
| All | 18 | +10.039 | -0.347 |
| All usable, same signature | 4 | +7.013 | -1.566 |
| All usable, multiple signatures | 1 | +31.688 | +0.000 |
| Some members unusable | 13 | +9.583 | +0.133 |

Initial HMM search remains on and downstream profile refinement is off. R is
conditional at P0/C0; C is conditional at P0/R0. These are not C-by-R interactions
or profile-refinement effects. The multiple-signature point is a single family,
not a replicated domain-order effect; no stratum-specific significance,
default change, causal explanation, FAS validation, independent confirmation
or superiority to full OrthoFinder is established.

## Execution And Checks

[Actual commands](swiss_ordered_pfam_execution_20261007_v1.json) retain ONE
invocation per stage using Python 3.10.13 -I -B/Biopython 1.87 and one-thread
library limits. Feature construction completed before its full independent
readback, which completed before prediction counts entered the projection.
The final exact-rational reader checks all 54 unchanged integer family count
rows, all 12 scores, eight contrasts, both TSVs and all 20 human numeric rows.
Failed R1 timing remains explicitly ineligible and unadmitted.

- Features: [c002bf9d, 331,517 bytes](swiss_ordered_pfam_features_20261007_v1.json).
- Feature readback: [275b8972, 8,669 bytes](swiss_ordered_pfam_feature_readback_20261007_v1.json).
- Projection: [f3389ad9, 40,391 bytes](native_swiss_ordered_pfam_20261007_v1/report.json).
- Rational/table readback: [db942129, 4,452 bytes](native_swiss_ordered_pfam_readback_20261007_v1.json).

Diagnostic postprocessing elapsed times/maximum RSS: features 0.25 s/41,472 KiB,
feature readback 0.22 s/41,472 KiB, projection 0.05 s/15,360 KiB, score readback
0.06 s/16,896 KiB; all exit 0 and zero swaps. These are shared-host observations,
not new admitted inference timings, isolated efficiency or comparable tool speed.
JSON/Biopython/csv parsing is shared. Annotation/raw-count verification is
inherited, not reproduced from original annotation files/reference trees here.

This closes a bounded ordered-annotation error-analysis gap, not literal full
domain architecture, fragment or ancestral-history truth. Families and earlier
outcomes are development-exposed; no new uncertainty or independent validation.
This supplement does not alter the frozen 21-page manuscript/review component
or establish whole-study publication readiness.

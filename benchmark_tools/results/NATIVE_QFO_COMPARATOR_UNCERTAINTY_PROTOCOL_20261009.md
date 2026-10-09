# Native-Cell SwissTrees Comparator Sensitivity Protocol

This analysis addresses the publication goal's paired method-difference
uncertainty requirement for the already admitted native QfO cells. Their point
estimates have been inspected. This is a retrospective, development-exposed
conditional sensitivity analysis, not prospectively untouched confirmation,
selection-adjusted inference, a new primary endpoint or method selection.
Freeze this protocol and tested implementation before calculating its intervals.
Do not tune a method or select cells after seeing the intervals.

## Verified Feasibility And Inputs

The [actual feasibility check](native_qfo_comparator_uncertainty_feasibility_20261009_v1.json)
checks all four native count audits against the two retained OrthoFinder modes:
exactly 18 common families, 563 represented proteins with disjoint family
membership, equal reference identity, and equal positive/negative truth totals
in each family. Recomputed macro precision/recall and harmonic F1 equal each
native count aggregate within 1e-14. This is direct retained-record checking,
not new raw scoring or native admission.

None of the four native family-record arrays equals an old eight-method array.
Consequently the original 24-endpoint comparator intervals cannot be relabeled
as these comparisons. The retained 42-endpoint factorial intervals describe
internal OrthoHMM contrasts, not OrthoFinder comparisons. Preserve both earlier
analyses and their correction inventories unchanged.

Use the exact six file identities in the feasibility receipt:

- Corrected eight-method comparator result: SHA256
  `a599749d66433ec211ed1ab0e3a6a2a75eb99cb0dc57abc47c2d267930cbbe34`.
- Four-cell native count/interval binding: SHA256
  `e8f554a178841da24a1892a7e42f8680f9fe8e4827dec6d6bfe4643cdb890bac`.
- Native P0/C0/R0 count audit: SHA256
  `4a49225b50de31f153c1b701f866db7acdfdbefa3d6398339ef5bf8711b2b561`.
- Recovered P0/C0/R1 count audit: SHA256
  `af4276cd759412a5a10ba04012a6f742f147c28611a6cff83adee463f85b0471`.
- Native P0/C1/R0 count audit: SHA256
  `7ed2ea1a7a94d55b28ecbadd7bd2be8e99df64b9fab3c266ebd987e8dbcfc479`.
- Native P1/C0/R1 count audit: SHA256
  `d5bafed830bf93e46b810728df14ffa34185d58ffe3eac78a693b0159d3d0a40`.

Require the current terminal-failure snapshot to have precisely these four
admitted cells and missing outcomes for indices 9/11/12. No partial successful
scorer output from a failed assessment is admitted. Compare every family
member list, reference identity and truth total, not just macro agreement.
Reuse original count conversion/statistic helpers without inventing another
pseudocount rule. Verify all stored family statistics and aggregate arithmetic.

## Fixed Comparisons And Missing Cells

Use eight factorial configurations in fixed P/C/R binary order. For each,
calculate candidate minus each of `orthofinder_3_1_5_full` and
`orthofinder_3_1_5_sequence_only`, for F1, PPV and TPR: 8*2*3=48 planned endpoints.
Only the four admitted native cells above have estimable contrasts. Preserve
all other endpoints as unavailable, not zero, and do not shrink the correction
family to the 24 currently estimable endpoints.

P1/C0/R0 is absent from this supplied fresh-cell snapshot; its retained
historical high-sensitivity benchmark is not missing globally and must not be
substituted into a new native cell. P0/C1/R1 failed before inference, P1/C1/R0
failed assessment, and P1/C1/R1 failed inference. Preserve those distinctions.
Initial sensitive HMM search remains on throughout. P toggles downstream profile
refinement; R also changes group-clique to native resolved-pair semantics.
Sequence-only OrthoFinder is the retained pre-phylogenetic MCL-checkpoint
diagnostic, not a separate timed final inference method.

## Fixed Calculation

Reuse the established family-bootstrap numerical kernel where compatible:
100000 multinomial draws of 18 families with replacement, PCG64 seed 20260920,
equal family probabilities and identical multiplicities across every available
native cell and comparator. Use the exact common canonical family order.
Regenerating the existing deterministic multiplicities is computation, not a
new independent dataset or permission to claim new independent confirmation.

Recompute weighted mean family PPV/TPR and their harmonic F1 in each replicate;
never average per-family F1 or treat gene pairs as independent. Report raw
0-to-1 differences, nominal 95% paired percentile intervals and intervals with
Bonferroni quantiles 0.05/96 and 1-0.05/96 across the fixed 48-endpoint family,
using NumPy linear interpolation. Report all 18 per-family differences and
descriptive wins/ties/losses with absolute tie tolerance 1e-10. Preserve every
contrast regardless of sign. No optional stopping, favorable seed choice or
rounding-based count matching.

Test inventory, pin, statistic and missing-endpoint refusals before production
execution. Independently recompute aggregate differences and inspect all
reported intervals/correction counts. Keep raw native serialized endpoints
distinct from full-precision count-reconstructed arithmetic. Use a fresh
output directory and commit/push prospective source/tests before executing.

## Interpretation And Boundaries

The 48-endpoint correction covers this explicitly separate retrospective
family, not all publication analyses or past tuning. Earlier 24/42-endpoint
analyses remain at their original scope. Eighteen reference families provide
limited conditional bootstrap resolution; disjoint proteins do not guarantee
biological exchangeability, exact percentile coverage or freedom from merged-
prediction correlations. Finite tail Monte Carlo resolution remains a limit.

An interval containing zero does not establish equivalence. An interval excluding
zero is conditional evidence for this exposed comparison, not generalization,
global superiority, an initial-HMM causal effect or a new recommended default.
No interval transfers to VGNC, TreeFam-A, GO/EC/FAS, the secondary mean, resource
timings or failed/missing native results. No native inference/scoring/retry,
public-source search, host/service change or additional timing run is authorized
by this statistical protocol. Full publication completion remains unproved.

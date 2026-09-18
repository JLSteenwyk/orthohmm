# Original-Release QfO Factorial: SwissTrees

## Scope And Validation

All eight original-release P/C/R cells have independently admitted scores.
P is profile expansion plus its downstream refinement, C is satellite_v2
candidate expansion, and R is the native inferred-pair strategy rather than
group-derived pairs. P-off still uses the initial HMM search. This is not an
HMM-free comparison or a new direct comparison with OrthoFinder.

The [frozen protocol](QFO_FACTORIAL_PROTOCOL_20260917.md) was applied without
retuning: 100,000 shared family-bootstrap draws, PCG64 seed 20260922, 18
SwissTrees families, 12 simple effects and two C-by-R interactions, each for
F1, precision and recall (42 adjusted endpoints). Each replicate recomputes
macro precision/recall and their harmonic mean, not the mean of family F1s.
Native raw counts use the audited half-count plus one prior.

Job 21725 completed successfully in 53 seconds on two shared-host CPUs.
A fresh execution of the frozen count auditor reproduced its entire report
exactly: eight cells share 10,765 reference relations with identical labels
and membership, and no genes overlap across the 18 reference families.
An independent calculation directly from TP/FP/FN counts, separately built
contrast weights, weighted family means and harmonic F1 reproduced all
point estimates, both interval sets, family differences and win/tie/loss
counts. Maximum absolute numerical discrepancy was 3.33e-16.

## Interpretation

- All 14 multiplicity-adjusted F1 intervals include zero. This neither
  establishes improvement nor equivalence.
- All four R contrasts increase precision and decrease recall, with adjusted
  intervals excluding zero in both directions. Observed precision gains are
  30.02-34.02 percentage points; recall changes are -4.56 to -6.81 points.
  Observed F1 gains are 10.40-13.12 points, but their adjusted intervals
  include zero. R changes prediction semantics and inference, not just group
  splitting, so the effect must not be attributed to splitting alone.
- With R off, C increases recall at both P levels; adjusted intervals exclude
  zero. With P off and R on, C has a small positive precision effect with an
  adjusted interval above zero. No C F1 benefit is established after adjustment.
- All 12 P-effect adjusted intervals include zero. This is not evidence that
  initial HMM search has no benefit: that search remains in every arm.
- Both C-by-R interactions have adjusted intervals including zero for all
  three metrics. Nominal evidence for the P-off F1 interaction does not
  survive the prespecified adjustment. Do not infer an interaction by
  comparing which simple effects are individually significant.

Only 18 development-exposed families are available. Shared evolutionary
history and merged predictions may violate family exchangeability; these
intervals do not adjust for development selection. They do not describe
uncertainty for other QfO challenges or the secondary six-metric mean.
The corrected-release reruns remain separate and unfinished; these results
must not be relabeled as corrected-release evidence or general superiority.

## Evidence

- [Eight admitted cells](qfo_factorial_admission_inventory_20260918.json),
  SHA-256 `4bf00bb5eba75335793c546a57daa0cbc378d8afd379b6dd306fd87006b019c1`.
- [Audited family counts](qfo_factorial_swiss_counts_20260918.json),
  SHA-256 `c513f9864255179bc1059ac8675aae713b4b393c9a94f50463249fc5eb1ed02c`.
- [Machine-readable statistics](qfo_factorial_swiss_bootstrap_20260918.json),
  SHA-256 `097eb4589d182521f3868316bf45a681657139ba7d6345abb2499ef2b6ae4751`.
- [Generated full 42-endpoint table](qfo_factorial_swiss_bootstrap_20260918.md),
  SHA-256 `4d79c6d80b0f368072a95591d51c8f4e1907cdb080c5bf2a7af390a53a9fa7af`.

The frozen analysis executor is
`benchmarks/work/publication_qfo_factorial_uncertainty_v1` at
`9082ccee291176b8883884d80e80ff4817053b86`. The
[batch workflow](qfo_factorial_uncertainty_batch_20260918.sh) performs the
count audit before bootstrap and refuses existing output paths. Preserve
the original artifacts when reproducing in a fresh output location.

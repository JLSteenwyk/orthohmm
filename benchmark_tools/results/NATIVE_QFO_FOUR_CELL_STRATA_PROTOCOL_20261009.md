# Four-Cell Native SwissTrees Fixed-Stratum Extension

## Question And Exposure

Extend the existing three-cell descriptive error analyses to the admitted
P1/C0/R1 cell, so that downstream-profile effects at C0/R1 can be reported
using the same fixed bins. This advances publication goal 4.3; it is not a
method change, new primary endpoint, independent validation or performance
improvement claim. Initial HMM search stays on throughout. Aggregate results
and the four changed profile pairs are already inspected; this extension is
retrospective and development-exposed. Commit/push this protocol before any
new profile subgroup calculation. Do not tune or select cutoffs after scores.

## Fixed Inputs

Use these exact retained artifacts:

| Input | Relative Result Path | SHA256 |
| --- | --- | --- |
| Three-cell fixed bins/counts | native_qfo_three_cell_strata_20261007_v1/report.json | 55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5 |
| Executed independent fixed-bin reader | native_qfo_three_cell_strata_readback_20261007_v2.json | 6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969 |
| Three-cell model-distance projection | swiss_model_divergence_strata_20261007_v1/report.json | ddfa27757c4c256df9c881079ac352c37d7e17646f75fbd417b624e7980be78f |
| Executed independent distance reader | swiss_model_divergence_strata_readback_20261007_v1.json | 186c0151137fd825513ef36c27df595b1a0f1cd4b2b4bf68de4b9665ad328d0e |
| Admitted profile-cell count audit | native10_allocated_swiss_counts_20261007_v1.json | d5bafed830bf93e46b810728df14ffa34185d58ffe3eac78a693b0159d3d0a40 |
| Executed independent profile count reader | native_qfo_profile_swiss_readback_20261007_v1.json | 3f62b487014cb9ba99ea36b5adc15535f62c7e39c9a99f3cb4eb8c537523005d |
| Four-cell admission/failure snapshot | native_qfo_terminal_failures_20261009_v1/report.json | 863ebc334dbdb5a41d1a9e9bba56be82896f0ff26c98366047928c7ef802966f |

Recheck direct artifact and source identities, including the new protocol and
prospective consumer. Reuse prior independently checked records; do not rerun
the old readers, scoring, HMM search, IQ-TREE, feature extraction or raw native
relation/proteome/tree/database parsing. Direct projection provenance is not
a renewed transitive native-output admission.

Require equal canonical memberships for all 18 families/563 represented genes
across the fixed-bin report, distance report and profile count audit. The old
three-cell family rows must agree across reports. All four cells must retain
the same positive/negative truth totals and 10,765 classified relations.
Check each cell's identity, original admission and raw-count binding against
the actual snapshot; recovered P0/C0/R1 timing stays failed/ineligible.
Validate all stored family prior statistics and full-precision macro scores.

## Frozen Bins And Contrasts

Reuse exactly the existing 11 sequence, five domain and four duplication bins,
plus the three model-distance bins: 23 suite/bin identities. Reuse every
membership list, cutoff and missing-feature designation; no new size proxy,
distance fit, bin cutoff, outcome ordering or fragment classification is
introduced. The model-distance lower/equal-median and higher bins stay at their
already frozen ordinary median of family medians. Repeated all-family and other
overlapping memberships remain audit rows, not independent findings.

Cell order: P0/C0/R0, P0/C0/R1, P0/C1/R0, P1/C0/R1. Contrast order:

1. R_at_P0_C0: P0/C0/R1 minus P0/C0/R0.
2. C_at_P0_R0: P0/C1/R0 minus P0/C0/R0.
3. P_at_C0_R1: P1/C0/R1 minus P0/C0/R1.

The new third contrast holds C0/R1 fixed; it does not identify profile effects
at other settings, candidate/reconciliation interaction or total HMM benefit.
R changes group-clique/native-pair semantics as well as reconciliation. Other
fresh unavailable configurations remain unavailable, not added from historical
methods or inferred from these contrasts.

## Statistic And Verification

Keep all integer TP/FP/FN/TN counts. Reuse the unchanged official halving and
unit-prior arithmetic: PPV=(TP+2)/(TP+FP+4), TPR=(TP+2)/(TP+FN+4) from raw counts.
Within each bin average family PPV/TPR and take harmonic F1, not pooled-pair
statistics or mean family F1. Empty bins retain empty status and NULL/NA scores
and differences, never zero or NaN.

Write a fresh separately named report and complete TSV/Markdown tables with
72 cell/family rows, 92 cell/bin rows and 69 contrast/bin rows. Reproduce all
69 prior cell/bin statistics and 46 prior contrast/bin statistics within
absolute 1e-12 at their actual original identities, not by schema renaming.
Keep human scores in percent and differences in percentage points; machine
records remain raw 0-to-1 values. Report every new profile bin regardless of
sign, including negative and neutral outcomes.

An independent standard-library reader must import no new exporter/projection
kernel. Use exact rational arithmetic on integer counts to verify all family,
bin, difference and human/TSV values plus inventories and input/output/source
identities. Tests cover prior preservation, mixed counts/memberships/admissions,
empty bins, priors, non-pooled aggregation, output tampering, scope inflation
and refused overwrite. Commit/push tested prospective sources before selected
generation at a fresh output. Preserve any failure; diagnose it before any
prospective recovery rather than silently overwriting or automatically retrying.

## Boundaries

No new bootstrap, subgroup interval, multiplicity-adjusted significance,
confidence-based bin selection, default promotion or native score/resource
admission occurs. Existing aggregate intervals do not become subgroup intervals.
Label-independent fixed features support descriptive associations only.
Length and missing explicit fragment flags do not establish fragment truth;
Pfam types/repeats are incomplete architecture descriptors; reference-node
duplication fractions are not ancestral duplication history. Estimated WAG+G4
branch-path distances include alignment/model, sampling and paralogy effects,
not dated evolutionary time or a known genealogy.

Auxiliary postprocessing durations are not OrthoHMM inference costs. Timings
remain shared-Threadripper observations with unknown, potentially tool-dependent
contention effects. No DGX, quiet window, service change or renewed routine
approval is required. Keep prior reports, figures, manuscripts and archives
unchanged; integrate the new descriptive findings later at their honest scope.
Continue the full seven-part publication goal afterward, not completion at this
projection milestone.

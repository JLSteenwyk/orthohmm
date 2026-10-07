# Native SwissTrees Three-Cell Fixed-Stratum Projection

## Question And Scope

Extend the completed two-cell descriptive error-stratum analysis to the newly
admitted candidate-expansion cell. This advances publication goal 4.3 without
altering OrthoHMM, choosing cutoffs from outcomes or rerunning a benchmark.
All three native cells retain the initial HMM search and have downstream
profile refinement off: p0_c0_r0 (6), p0_c0_r1 (7), p0_c1_r0 (8).
Reconciliation and candidate expansion are separate conditional contrasts
against cell 6. These three cells cannot identify their interaction.

This protocol is committed and pushed before calculating any candidate
subgroup scores. Existing aggregate candidate results have already been
inspected; this is development-exposed descriptive analysis, not a prospectively
unseen confirmation. The sequence, domain and duplication bins were frozen
before this extension and are not changed or ranked to select favorable bins.

## Inputs And Frozen Memberships

Use exactly these retained artifacts and their bound original sources:

- Sequence report4650f614 and independent readback8d0356bc.
- Domain reportd3235c71 and independent readback9e77b67a.
- Duplication reporte0ecbc9d and independent readback84eb3d9f.
- Sequence feature/membership artifact91236327.
- Candidate count audit7ed2ea1a and independent rational readback07941ba9.

Recheck current bytes of these direct inputs and the original sources named by
their retained readbacks; do not rerun those readers. Reuse memberships from
the verified reports. Require exact agreement for all 18 families and 563
canonical accessions, including the duplication report's already validated
alias-to-entry resolution. Do not reread native relations, FASTAs, annotation
JSONs, mapping or database files, original trees, candidate merges or old
bootstrap draws. This projection's checks are direct artifact/source checks,
not renewed raw evidence validation or full transitive provenance admission.

Preserve the 11 sequence bins, 5 domain bins and 4 duplication bins, with
each report's original membership lists. The all-family bin is repeated in
each suite for auditability, not counted as three independent findings.
Require bin membership consistency across the old two cells, agreement with
explicit domain/duplication bin dictionaries, and all old family counts and
statistics identical across the three reports. Candidate represented-gene
sets must exactly match each reference family; family total classified pairs
and TP+FN reference-positive totals must remain identical in all cells.

Sequence descriptors are length/composition proxies, not calibrated divergence
or literal fragment truth. Domain bins describe Pfam annotations, not complete
validated architectures. Duplication bins are reference-derived informative-node
fractions, not reconstructed ancestral duplication histories; retain the
original default-S limitation. Absence of fragment labels is not completeness.

## Statistics And Outputs

For each cell/family retain original integer TP, FP, FN and TN counts without
prior. Apply the QfO convention: TP/2+1, FP/2+1, FN/2+1. Calculate family
precision and recall, then macro-average precision/recall over the bin and
take their harmonic F1. Do not pool pairs, average family F1, remove TN from
the preserved counts or treat dependent pairs as independent replicates.

Publish all 54 cell/family rows, all 60 cell/bin rows and all 40 contrast/bin
rows (20 bins times two conditional contrasts). Empty bins have explicit
empty status and NULL/NA metrics and differences, not zero, NaN or omission.
Rows retain unit-interval values; human tables display percentages and
contrasts in percentage points. Verify reproduction of every prior cell/bin
statistic within absolute 1e-12; preserve all negative and neutral effects.

Write a new versioned report, complete TSV tables and a human-readable table;
never overwrite the two-cell reports or their bound exporters. An independently
implemented stdlib reader imports no exporter or prior scoring/projection
kernel. It rechecks direct input/source/output identities, family/member/bin
inventories, the complete TSV/JSON tables and human table values. It recomputes
all statistics from integer counts using exact rational arithmetic, with
absolute 1e-12 comparison to exported numeric values. Include focused fixtures
for count conversion, non-pooled aggregation, empty bins, membership/count/bin
tampering, scope inflation, output alteration and refused overwrite.

No new bootstrap draws, attached subgroup intervals, multiplicity-adjusted
significance, cutoff optimization, score/timing admission, confidence or causal
claim. The recovered cell 7's failed timing remains ineligible. These scripts'
postprocessing durations are not inference resources. The original timings
remain shared-host observations with unknown, potentially tool-dependent
contention effects; do not imply isolated efficiency or repair failed timings.

## Execution And Continuation

Commit and push protocol, then tested new sources before selected projection
into fresh native_qfo_three_cell_strata_20261007_v1. Use retained scientific
Python 3.10 in isolated stdlib mode and one fresh independent readback. Retain
any actual execution failure without silently overwriting or automatic retry.
Commit/push resulting evidence and focused result summary after validation.
Existing native job23902 and review23910 continue unchanged; do not resubmit,
restart or release them again. Continue the full publication goal afterward,
not completion at this descriptive milestone.

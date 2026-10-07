# Prospective Ordered-Pfam Annotation Error Analysis

This is a bounded continuation of original publication goal 4.3, not a new
accuracy campaign, method-selection experiment or claim of complete biological
architecture. Existing SwissTrees domain strata use type/repeat counts; this
analysis uses retained domain order. Freeze this protocol and tested sources
before constructing selected descriptors or projecting prediction counts.

## Frozen Population And Inputs

Use all 18 corrected SwissTrees reference families and all 563 canonical members.
Reuse, without annotation extraction, realignment or filtering:

- `swiss_domain_annotation_inventory_20260917.json`, SHA-256
  `d5e269158c2fb603acc0805a1c140a88758342d7b8e32b3094750ade22fbc06c`.
- `native_qfo_swiss_domain_strata_20261006_v1/report.json`, SHA-256
  `d3235c712e711f005b3182a3621c0079f8a85c3c7abd5db016d5e8ea95884cd6`.
- Its independent annotation/readback receipt, SHA-256
  `9e77b67ad97e6824e9442691be7a05ba1552d62ffe696ced6034a996ede1945d`.
- Retained model-distance feature report, SHA-256
  `ceb4ca3675836e1815fcc1f14c29c8a50947633a8a5e8dd6c67477bdd853065c`;
  its 18 admitted alignment bindings supply sequence IDs and ungapped lengths,
  not distances, outcomes or new phylogenetic labels.

The prediction-independent feature stage must not load the native count report.
Prior annotation verification is inherited, not rerun over all 78 annotation
files. Confirm all retained feature identities and canonical member sets.
Validate supplied Pfam coordinates/type/repeat counts; require annotation length
to equal the corresponding ungapped retained protein length. A missing member,
length mismatch or invalid instance fails the complete stage. No successful-
only subset supplies bins or scores. Every selected family must be retained.

## Descriptor And Bins

For every protein, sort all retained Pfam instances by `(start, end, domain)`.
Preserve repeat multiplicity. A usable ordered signature requires at least one
instance, `start < end` for every instance and strictly separated intervals:
each next start must be greater than the preceding end. Ties, touching,
overlaps/nesting and zero-span records are conservatively order-ambiguous.
This rule does not assume inclusive versus half-open coordinate widths.
Zero Pfam hits are a separate unusable annotation state, not biological absence.
Never resolve ambiguity by a domain-name tie-break, discard a competing domain,
impute an absent type or treat an unflagged annotation as complete architecture.

Retain per-protein status, ordered signature when usable, the complete ordered
instances and sorted type/multiplicity multiset. Per-family descriptors include
member/usable/ambiguous/zero-hit counts and distinct usable signatures.
The four primary reporting bins are fixed, with no learned cutoff:

1. `all`: every family.
2. `all_members_usable_same_signature`: every member usable and one signature.
3. `all_members_usable_multiple_signatures`: every member usable and at least two.
4. `some_members_unusable`: at least one member unusable, irrespective of the
   remaining members' signatures.

The last three bins partition all 18 families. These labels describe the
retained annotations, not experimentally verified completeness or domain loss.
Retain empty bins as NA; never redefine bins because they are empty or show an
unfavorable effect. No alternative cutoffs, percentiles or subgroup selection.

As an unscored feature diagnostic, count all within-family unordered pairs of
usable proteins with identical Pfam type/multiplicity multisets, and count which
have different ordered signatures. Retain zero-comparison cases explicitly.
This distinguishes observed order variation from different type content or
repeat count; it is not a separate score stratum, interaction or hypothesis test.
Within-species and paralog comparisons remain in this descriptor.

## Separate Outcome Projection

Only after complete feature construction and a separate readback may the
three-cell counts enter. Use the already admitted full count report
`native_qfo_three_cell_strata_20261007_v1/report.json`, SHA-256
`55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5`,
and readback SHA-256
`6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969`.
Retain all 54 original family count rows (TP/FP/FN/TN), member sets and
cell semantics: P0/C0/R0, P0/C0/R1, P0/C1/R0. Initial HMM search is on,
downstream profile refinement off. Preserve failed R1 timing as ineligible.

For each nonempty bin/cell, apply the original TP/2+1, FP/2+1, FN/2+1 prior,
mean family precision and recall, then harmonic F1. Do not pool counts or average
family F1. Retain all 12 score rows and eight conditional differences:
R at P0/C0 and C at P0/R0, in percentage points. Empty-bin statistics are null.
An exact-rational readback checks every displayed metric against retained
integer counts; a separately implemented feature readback checks all protein
states and bins, using all-pair interval-overlap checks rather than the primary
adjacent-interval rule. Shared JSON/Biopython parsing remains explicit.

## Scope And Failure Rules

No raw scorer, original tree traversal, bootstrap, HMM search, new inference,
annotation generation or OrthoFinder run. No new intervals/significance,
default promotion, outcome-selected threshold, interaction, calibrated
architecture, causal biology, independent confirmation or FAS validation.
The families and earlier outcomes are development-exposed. Ordered annotations
can reflect divergence, length, missing hits and annotation assumptions.
This advances an ordered-annotation error analysis; literal fragment, full
architecture and true ancestral-history evidence remain limited.

Retain all failures and original bytes. Use distinct fresh output namespaces,
commit/push tested sources before selected execution, and never automatically
retry failed analyses. This small postprocessing task does not gate native
23902/23910, require a quiet host or alter their sources/scheduler state.

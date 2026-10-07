# Profile Refinement: Where Four Predictions Changed

This chronological supplement explains every changed assessed SwissTrees
pair in the P1C0R1 minus P0C0R1 profile-refinement contrast. It supplements the
[conditional uncertainty result](NATIVE_QFO_PROFILE_UNCERTAINTY_RESULT_20261007.md);
old manuscript, figures and archives remain frozen. Initial HMM search is on
in both cells; this is not a total-HMM control or an OrthoFinder comparison.

## Observed Changes

The [selected trace](native_qfo_profile_localization_20261007_v1.json), SHA256
`0d09567de8828b46312aa280f0bdabe6bbfbfaa01f591cbbca9190241bdb5ae0`,
reconstructs the complete pre-reconciliation candidate partitions from the
admitted root-HOG files, checking original input-cluster hashes and input
universe. It checks observed node-table topology/event rules and streams
native pair files to verify actual inclusion. It finds four changed scored
relations, with no additions. All four endpoints were in a common candidate
family with a speciation pair-event before refinement; each pair spans two
different candidate families afterward. Thus these removals already follow
from candidate separation, not a new duplication exclusion at a shared LCA.

| Reference Family | Protein Pair | Change | Before Candidate | After Candidates | Software Localization |
| --- | --- | --- | --- | --- | --- |
| CASP | Q0DHM7 / Q8L7R5 | FP to TN | Family0061566 (11 genes) | Family0061640 (10) / Family0061018 (16) | Candidate separation |
| CASP | Q5N794 / Q8L7R5 | FP to TN | Family0061566 (11 genes) | Family0061640 (10) / Family0061018 (16) | Candidate separation |
| CASP | Q6K478 / Q8L7R5 | FP to TN | Family0061566 (11 genes) | Family0061640 (10) / Family0061018 (16) | Candidate separation |
| GH14 | Q10LG9 / Q8VYW2 | TP to FN | Family0115627 (7 genes) | Family0135675 (6) / Family0083578 (20) | Candidate separation |

The before-state CASP pairs share node G00009; GH14 uses G00005. Both annotated
pair-events are speciation. The after states have no within-family LCA for
these pairs; inventing an after-state duplication node would be incorrect.
All endpoint candidate-member hashes change. Three removed false positives
improve CASP precision, but the lost true positive worsens GH14. The full
18-family macro F1 difference is -0.320215 percentage points, with a
42-endpoint-adjusted interval [-2.112886, 0.513227]. Sixteen families tie.
Removing more false positives than true positives does not ensure an improved
macro statistic: the family-level changes receive equal family weight and
affect precision and recall differently.

## Independent Saved-Tree Check

Source was committed/pushed at d27ca3ae before one selected
[Newick readback](native_qfo_profile_newick_readback_20261007_v1.json), SHA256
`289ac094115f493b8607073b1395ca2d97493ba4a00d5f1a083ceb3aa030507f`.
This uses Biopython1.87/Newick LCAs and a separate streamed CSV root-membership
parser, not the primary node-table LCA code. It checks six saved family
checkpoints and 70 tree leaves, raw/rooted/annotated digests, identical
rooted/annotated topology, complete candidate members/hashes/sizes, all four
pair states and both distinct before-state LCAs. All four localizations are
candidate separation; no changed pair is a positive-paralogy exclusion.
Tree/checkpoint records are newly observed and digest-consistent, not
retroactively relabelled as originally inventoried admission evidence.

Trace source was committed/pushed at247033cd before its single execution;
exit0/38.50s/250248KiB RSS/zero swaps are retained in its
[execution receipt](native10_profile_trace_execution_20261007_v1.json).
The independent reader exits0/5.80s/36864KiB RSS/zero swaps, retained in its
[execution receipt](native10_profile_newick_execution_20261007_v1.json).
These are read-only diagnostic costs, not inference or isolated efficiency.
All59 focused invented-tree/corruption/helper tests pass in1.60s. Selected
sources/results are now bound and must not be edited or automatically rerun.

## Interpretation And Limits

The frozen pipeline expands profiles, combines those edges with RBH-derived
edges, reclusters, recomputes final singleton edges and clusters again before
phylogeny. Membership need not grow monotonically when extra evidence changes
graph clustering. The observed separation is consistent with that workflow;
these artifacts do not identify the individual profile hit, graph edge or
clustering decision that caused it. No earlier intermediate search/profile
history has been reconstructed here.

The inferred species-tree byte hashes differ between the two cells. Separation
before reconciliation explains these four software absences independently of
an after-state LCA, but does not isolate the overall effect of species-tree
changes. Saved annotated trees do not establish true duplication histories,
correct rooting, calibrated confidence or biological correctness. Checks use
the same development-exposed outputs and are not independent biological
confirmation. The recovered reference timing remains failed/null/ineligible.
There is no new inference, alignment, rooting, reconciliation, conversion,
scoring, FAS sample or bootstrap, and no default-promotion, general superiority
or publication-readiness claim. The native factorial remains incomplete.

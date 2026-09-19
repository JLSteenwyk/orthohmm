# Corrected QfO Parameter Neighborhood

## Scope And Exposure

Execute the same six one-at-a-time variants fixed in
`PARAMETER_NEIGHBORHOOD_PROTOCOL_20260916.md`, whose SHA-256 is
`6f6f98eb2053348499dff7656d90d5ce4f38127b0de966bdb2cdeedfde17eb79`.
The original protocol deferred QfO until reproducible baseline admission.
That prerequisite is now met by the corrected replay and complete factorial.
This supplement fixes QfO-specific provenance and uncertainty before any
neighborhood variants are run or scored. The baseline and factorial outcomes
are already known; this is development-exposed sensitivity analysis, not
independent confirmation, prospective method selection or a new default search.

## Baseline And Fixed Panel

Use the corrected 984,137-protein, 78-species input and frozen CPU method.
The full-pipeline baseline is corrected factorial `p1_c1_r1`: profile
refinement, satellite_v2 candidate expansion and inferred phylogeny enabled.
`qfo_parameter_neighborhood_plan_20260919.json` binds the existing replay,
candidate manifest, native admission and score admission. Verify those
records and their transitive runtime/input/checkpoint records before use.
The checked replay's profile-refined partition is equal to the corrected
native high-sensitivity partition; this is not an assumption of equality
between different input releases.

Keep the original arm order: control, cpm_low, cpm_high, norm_low, norm_high,
margin_low, margin_high. Resolutions are 0.08 and 0.12 versus 0.1; min_norm
is 0.024 and 0.036 versus 0.03; min_margin is 1.2 and 1.8 versus 1.5.
No joint perturbations, added seeds, adaptive ranges, outcome-selected
repeats or replacement of unfavorable variants. All other settings remain
fixed. No use of YGOB outcomes to choose settings.

## Execution And Admission

Candidate-threshold arms reuse the admitted corrected profile-refined seed
partition and numeric hit checkpoint. First repeat unchanged candidate
expansion and require exact baseline candidate-partition and merge-trace
bytes. Use the existing scoped override of the candidate engine, recording
both applied parameters and the unchanged nominal wrapper report. Restore
the engine on failure. Do not edit the scientific package.

CPM arms must rerun graph clustering, singleton handling, profiles and
refinement from the fixed initial-hit checkpoint. First require the unchanged
CPM control to reproduce all four corrected replay stage partitions; do not
reuse the baseline profile-refined partition for changed CPM resolution.

For every changed candidate partition, infer its species tree and rerun
rooting, reconciliation and membership constraints. Reuse raw alignments or
gene trees only under the existing exact membership, sequence and tool/hash
checkpoint rules. A fixed supplied tree is not a substitute. Where candidate
partitions, constraints, inputs and complete commands are identical, record
the equivalence explicitly before sharing native outputs; never infer such
equivalence from aggregate group counts or scores.

Run through the scheduler on the workstation, not the occupied DGX. Candidate
preparation uses 2 CPUs; grouping/profile and reconciliation use the existing
32-CPU limits and frozen numerical-thread controls. Capture actual allocations,
commands, exit codes, checkpoints and incremental time/memory. These cached,
shared-host observations are not controlled end-to-end speed comparisons.
Do not stop unrelated jobs or change the active primary comparator chain.

Independently validate native outputs and cross-species pair conversion before
official scoring. Retain every failure and its stage; a missing or pending
run is neither a failure nor a zero. No partial output is admitted as complete.
The protocol/plan freeze is not authorization to bypass unimplemented runner
or admission checks. Freeze and test those executors before submission.

## Endpoints And Uncertainty

Report all six official QfO endpoint summaries and their native axes, prediction
counts, input/reference coverage and conversion losses. The unweighted
six-summary mean is project-defined and secondary. Do not reuse historical-
release scores or conflate pair counts with covered proteins.

For SwissTrees use the same 18 reference families and independently audited
native counts. Use 100,000 shared PCG64 family-bootstrap draws, seed 20260925,
and linear percentile intervals. Within each draw recompute macro precision
and recall, then their harmonic-mean F1; do not average family F1. Compare
each variant minus the unchanged full-pipeline baseline and report family
wins/ties/losses. Report nominal 95% intervals and Bonferroni intervals over
all 18 planned endpoints (six variants times F1/precision/recall), retaining
that denominator if any variants fail. Failed variants have no intervals;
if the control is unavailable, no paired contrasts are estimable.

Do not treat dependent protein pairs as independent resampling units. No
TreeFam family intervals or uncertainty on the heterogeneous secondary mean
are invented. These intervals remain conditional on development-exposed
families and are not adjusted for earlier method selection. An interval
including zero does not demonstrate equivalence or general robustness.

Retain every completed, failed and pending arm in the final inventory. No
variant is promoted as a new default and no competitor superiority is inferred
from a favorable neighborhood result. This panel addresses local parameter
sensitivity, not the remaining matched-search or controlled-scaling gaps.

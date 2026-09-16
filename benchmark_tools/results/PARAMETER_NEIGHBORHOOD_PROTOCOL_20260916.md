# Limited Parameter Neighborhood

## Purpose And Scope

Post-development sensitivity analysis, specified before running the neighborhood
or viewing its scores. Existing QfO/OrthoBench results, ablations and traces have
influenced development; this is not independent confirmation. Freeze a bounded
one-at-a-time neighborhood to describe stability of the publication method,
not hill-climb or select the best variant. No new default will be promoted from
this panel. Any later method change requires fresh independent confirmation,
not retuning against the already inspected YGOB results.

## Fixed Panel

All other frozen HMM search, profile, refinement, candidate and phylogeny
settings stay unchanged. Each variant changes exactly one scalar by20%:

| Arm | Parameter | Baseline | Variant |
| --- | --- | ---: | ---: |
| cpm_low | CPM resolution | 0.1 | 0.08 |
| cpm_high | CPM resolution | 0.1 | 0.12 |
| norm_low | Candidate min_norm | 0.03 | 0.024 |
| norm_high | Candidate min_norm | 0.03 | 0.036 |
| margin_low | Candidate min_margin | 1.5 | 1.2 |
| margin_high | Candidate min_margin | 1.5 | 1.8 |

Include the unchanged control. No joint perturbations, adaptive grid expansion,
outcome-selected repeats or seed selection. All variants retain the HMM-centered
method, satellite_v2 constraints and full inferred-tree pipeline. Changes in
candidate families can change species-tree inference; that downstream effect
belongs to the end-to-end parameter contrast. Do not substitute a fixed supplied
tree without explicitly labeling a separate diagnostic.

## Execution And Admission

OrthoBench runs first. Reuse the verified normalized-hit cache. For candidate
threshold variants, reuse the byte-admitted profile-refined seed partition.
For CPM variants, repeat the graph/profile/refinement stages at the changed
resolution; do not reuse the baseline profile-refined partition. Reuse raw gene
trees only under the existing exact-membership/input-hash checkpoint rules.
Infer the species tree and rerun rooting/reconciliation/constraints for every
changed candidate partition. Preserve cache/source/tool/input/output hashes,
actual parameters, failures, coverage and incremental resource measurements.

Candidate preparation first requires the unchanged control to reproduce both
the saved candidate-partition bytes and the complete merge-trace bytes. The
four threshold variants use the frozen engine with an analysis-only scoped
call override. Record the wrapper's unchanged nominal profile report separately
from the actual applied parameters and verify exactly one engine invocation.
Restore the engine after each call, including failure. No production source or
default is edited. This preparation does not constitute completed reconciliation
or accuracy evaluation, and it does not execute the two CPM variants.

Require complete native group integrity and correct output semantics before
scoring. Failures remain failures, not zeros or exclusions hidden from tables.
Preserve expensive outputs and avoid needless search/alignment reruns. Shared-node
preparation/reconciliation times are incremental observations, not controlled
end-to-end efficiency evidence.

QfO receives the same six-variant panel only after its unresolved reproducible
baseline is established. Do not interpret the current affinity experiment as
proof that a preferred historical partition should be selected.

## Endpoints And Reporting

For OrthoBench compare each full-pipeline variant with the unchanged full-pipeline
OrthoHMM control. Report F1, precision and recall, complete coverage and family
wins/ties/losses. Use20,000 paired RefOG bootstrap draws, seed20260918, retaining
the actual full-reference weighted statistic and low-certainty conventions.
Bonferroni adjustment includes all18 planned endpoints (six variants by three
metrics), even if a run fails. Report missing intervals explicitly; no reduced
multiplicity count after observing results. Do not treat gene pairs as independent.
Candidate-only and intermediate-stage scores, if shown, are descriptive diagnostics,
not additional selected efficacy tests.

Show individual QfO endpoints and coverage, preserving each benchmark's native
semantics and applicable uncertainty procedures. The project-defined six-metric
mean is secondary; do not invent a common independent sampling unit across
heterogeneous QfO endpoints. State unresolved uncertainty limitations explicitly.

Retain all neutral/adverse outcomes. A confidence interval including zero does
not prove equivalence or robustness, and a positive difference on this exposed
panel does not establish superiority over OrthoFinder or generalization.

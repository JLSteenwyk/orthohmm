# Three Native Collector Paths Validated

The remaining two private-controller collector-v5 fixtures completed sequentially.
The [OrthoHMM phylogeny submission](threadripper_private_collector_submission_22374.json)
ran first; the [OrthoFinder submission](threadripper_private_collector_submission_22375.json)
had an explicit successful-completion dependency on it. Neither identity was
retried. These are four-species/16-gene installation fixtures, not production
timing runs or new accuracy results.

| Native path | Job | Terminal allocation/batch/step | Groups | Native phylogenetic pairs | Replay and retained prediction comparison |
|---|---:|---|---:|---:|---|
| OrthoHMM high sensitivity | 22373 | COMPLETED 0:0 | 3 | Not applicable | Pass (previous turn) |
| OrthoHMM satellite_v2, inferred phylogeny | 22374 | COMPLETED 0:0 | 3 | 36 | Pass |
| Full OrthoFinder 3.1.5 | 22375 | COMPLETED 0:0 | 3 (pre-phylogeny checkpoint) | 36 | Pass |

The [phylogenetic OrthoHMM audit](threadripper_private_collector_22374.json) and
[full OrthoFinder audit](threadripper_private_collector_22375.json) bind local
raw replay, native output validation, terminal scheduler observations and
canonical prediction comparisons against the retained previous fixtures.
Both before/after runtime checks pass in both jobs. OrthoHMM's three root
groups and native pairs match; OrthoFinder's checkpoint membership and native
pairs match. Input/species mappings match as well. Group labels and row order
are normalized by the existing fingerprint implementation; this comparison
does not claim every intermediate file is byte-identical.

Independent outcome auditing replays the raw collector evidence and validates
successful native formats separately. Both completion receipts report
`anchor_only_at_boundaries`. OrthoFinder's native species/sequence mapping
checks pass. Forty-one focused tests pass across replay, outcome auditing and
canonical fingerprints. No scientific configuration, dependency or helper
source changed during these jobs. The Slurm queue was empty after completion.

## Remaining Production Requirements

Together with [job 22373](THREADRIPPER_PRIVATE_COLLECTOR_20260928.md), this closes
the small-fixture integration check for all three native collector paths under
the private deployment. It does not establish full-scale correctness,
whole-run quiet-host eligibility, complete controller branch coverage or
final job-teardown resource accounting. The shared-host fixture times are
retained only as diagnostics and must not enter the 27-run timing panel.

Next, resolve the outstanding production resource-accounting and environmental
admission requirements, including a defensible quiet window. Continue other
publication analyses if the host cannot provide one. No production identity
has launched, no unrelated workload was stopped, and no DGX access occurred.

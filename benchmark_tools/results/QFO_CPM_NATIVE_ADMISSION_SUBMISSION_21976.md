# CPM Native Phylogeny Admission

## Status

Submitted September 19, 2026 at 10:27 EDT as serial array `21976_0..1`.
Both tasks were confirmed PENDING with dependency reasons immediately after
submission. No CPM native result, accuracy score or timing comparison is
admitted by this submission record.

Executor commit: `cdbf61402d7d77512bc86005e059c18a353dcc62`, pushed to
`origin/main`. Frozen worktree:
`benchmarks/work/publication_qfo_cpm_phylogeny_admission_v1`.

| File | SHA-256 |
| --- | --- |
| `benchmark_tools/admit_qfo_cpm_phylogeny.py` | `989248b766677494d5377f8c3b3e0d4d1642b18e5ec13e97cde53294dd3277b0` |
| `benchmark_tools/results/qfo_cpm_phylogeny_admission_batch_20260919.sh` | `5626027707ef3a651a00b06441cf2d5f1a4cc7b10284387aea73788bf479b156` |

## Submission

```bash
sbatch --parsable benchmark_tools/results/qfo_cpm_phylogeny_admission_batch_20260919.sh \
  /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/publication_qfo_cpm_phylogeny_admission_v1 \
  cdbf61402d7d77512bc86005e059c18a353dcc62
```

The scheduler confirms two serial tasks, each requesting 2 CPUs, 64 GiB,
four hours on `bizon`, without requeue. Dependencies are both
`afterany:21972` and `aftercorr:21972`: wait for the full producer array to
be terminal, then require successful completion of the corresponding task.
This avoids comparing changing whole-array accounting snapshots. Index 0
is `cpm_low` (resolution 0.08); index 1 is `cpm_high` (0.12).

Reports will be written exclusively to
`benchmarks/work/qfo_cpm_native_admission_21976_<index>.json` and logs to
`benchmarks/work/qfo_cpm_native_admit_21976_<index>.log`.

## Validation Scope

The validator requires a completed 32-CPU producer task before reading its
outputs. It verifies the frozen producer commit/source, complete Python
helper inventory, exact execution provenance, tool paths, input hashes,
candidate-admission evidence and fresh candidate-admission report. The
full CPM-specific seed/candidate/constraint arm replaces the baseline arm
in the native validator's manifest; no baseline seed is silently reused.

Shared validators check successful process/artifact records, inferred-tree
provenance and taxon coverage, reconciliation rules and membership counts,
complete native group coverage, and cross-species native ortholog pairs
within their candidate families. Inputs, executable sources and native
artifact inventories are rechecked before admission. This does not
independently reconstruct gene trees or orthology events.

Successful status is `cpm_native_pairs_verified_unscored` with
`accuracy_evaluated`, `scoring_admitted`, and `publication_ready` all false.
Mapping, conversion, scoring and score admission remain separate steps.
Shared-host inference time is not admitted as controlled comparative timing.

## Verification

230 tests passed in 12.75 seconds across the new validator, CPM producer,
candidate admission, parameter native admission, native metadata/partition
validation, process validation, factorial native pair checks and launcher
checks. Tests include complete mocked admission orchestration and reject
wrong scheduler state/resources, changed source/inputs, incomplete helper
provenance, mismatched fresh admission, wrong CPM partitions, empty pair
output, and altered postflight evidence. These tests are not a substitute
for validating the eventual dataset outputs. Batch `bash -n` and staged
`git diff --check` passed.

## Progress Ledger

- Implemented, tested, committed and pushed native CPM admission; queued
  behind inferred phylogeny array `21972`.
- At 10:27 EDT, OrthoFinder score admission `21736_0` was RUNNING, with
  `21736_1` pending. Newly completed scoring is not yet an admitted result.
- Legacy BLAST `21713`, DGX timing task `21920_16`, and host-side recorder
  `21922` were RUNNING; timing task 17 remained pending. No DGX SSH/native
  output inspection was performed during its timing panel.
- Original TreeFam family trees/mapping remain unavailable; no additional
  family-level TreeFam uncertainty is enabled.
- Next: CPM conversion, QfO assessment, score admission, and integration
  into prespecified parameter uncertainty; retain comparator admissions
  and audit the DGX panel only after all timing tasks are terminal.
- Publication readiness remains unproven. No scientific defaults or
  endpoints changed, and unrelated working-tree changes were preserved.

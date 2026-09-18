# Corrected QfO Reconciliation Admission

The corrected reconciliation executor is frozen at
`2bf70fb27cc63edc7c49d38d3c7f7d09ddccca9b` in
`benchmarks/work/publication_qfo_corrected_reconcile_v1`.
This freezes the existing runner; it does not change the scientific core,
replay launcher, candidate algorithm or reconciliation settings.

`admit_qfo_corrected_factorial_cell.py` admits one successful reconciliation
array element using its **raw scheduler job ID** and index 0-3. No corrected
reconciliation has yet completed; this is implementation evidence only.

## Gates

- Require a terminal successful 32-CPU `bizon` job before reading candidate
  or reconciliation output. Reject running, failed and cancelled tasks.
- Independently reverify the admitted corrected candidate report, its
  scheduler identity, frozen preparation, all bound inputs and runtime.
- Require the exact frozen reconciliation executor and clean tracked source.
- Reconstruct and compare the complete execution provenance: relocated
  command, candidate admission, tool resolution, input records, scheduler
  task, working directory and source equivalence. Require corrected postflight.
- Verify the complete inference artifact inventory and hashes, then reuse
  the native metadata, provenance, tree, membership and RootHOG validators.
- Require the full 984,137-gene corrected universe and 78 distinct taxa.
- Validate native pair ordering, uniqueness, species ownership, candidate
  containment and agreement with the completion summary. Do not manufacture
  RootHOG clique pairs in place of native phylogenetic predictions.
- Recheck output inventory, candidate admission, source records and runtime
  before returning an admission. CLI writes only to a fresh result path.

Successful status is `corrected_qfo_native_pair_output_verified`.
`accuracy_evaluated`, `scoring_admitted` and `publication_ready` remain false.
This checks native consistency, not independent phylogenetic truth or
event-to-pair reconstruction. Corrected pair conversion and scoring adapters
are still required; historical results cannot fill those roles.

## Execution

After candidate admission, use the existing corrected reconciliation batch
with the frozen executor above. Once each reconciliation element is terminal,
submit `qfo_corrected_factorial_admit_batch_20260918.sh` with:

```text
ADMISSION_EXECUTOR ADMISSION_COMMIT INDEX RAW_RECONCILIATION_JOB
CANDIDATE_ADMISSION_PATH CANDIDATE_ADMISSION_SHA CANDIDATE_ADMISSION_JOB
FRESH_OUTPUT_PATH
```

The admission executor must contain the new checker and be frozen separately
after validation. The batch requests 2 CPUs, 64 GiB and four hours with no
requeue. Do not substitute the parent array ID for the raw element job ID.

## Validation Scope

Focused tests cover terminal/resource/task checks, complete provenance and
input identity, the full corrected identifier universe, duplicate members,
and the reused native-pair and runner gates. A live negative check against
running job21707 rejected admission before nonexistent candidate files were
accessed. The batch passes `bash -n`. A positive full-scale corrected native
admission remains untested until the actual upstream jobs finish.
The complete unit suite passed with 3,123 tests and one optional skipped
test in 63.05 seconds; this does not replace that pending native-data check.

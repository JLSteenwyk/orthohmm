# Explicit Helper-Runtime Recovery Admission Amendment

## Scope

Preserve the [original recovery protocol](QFO_CPM_CHECKPOINT_RECOVERY_PROTOCOL_20260923.md),
failed job 22081_1, completed continuation 22154, failed admission 22155 and
allocator failure 22158. Their reports/statuses/executors remain unchanged.
The [separate helper-complete control](QFO_CPM_HELPER_REFINEMENT_RESULT_20260930.md)
completed the original unchanged refinement-only runner in a reconstructed
runtime, producing the same retained partition and metadata except destination.
It did not itself admit a seed. This new read-only admission names that runtime
amendment rather than impersonating the failed original admission.

Only independent refinement execution uses this explicitly different helper
runtime. The completed original optimizer retains its actual historical runtime,
arguments, constructor adapter, input graph and native before/after observations.
The refined seed bytes remain the originally recovered bytes. No new parameter,
endpoint, output choice, optimizer/refinement attempt or accuracy evaluation.
No claim that the crash is repaired, its cause known or all runtimes equivalent.

Write a distinct admission status/output root. Original admission/candidate
validators keep their original meanings; do not change their default behavior
to accept a failed job. Seed admission alone must not release old dependencies.
An explicit recovery-aware candidate handoff is a separate next gate. Subsequent
candidate science/runtime must remain matched to the other parameter arms, not
silently switch to the helper runtime because it passed this refinement control.

## Required Evidence

Commit/push `admit_helper_cpm_recovery.py`, tests and this protocol before one
actual audit. Require this protocol's explicit SHA256. Fixed input digests bind
the original parent/preflight/failed admission and the two independently checked
new-runtime preparation/refinement receipts. Reject overwrite, drift, absent
proof or conflicting identities; retain failed preflight/audit attempts without
automatic retry. No native scientific execution is permitted by this audit.

1. Current Slurm accounting must show 22154 COMPLETED 0:0 on one CPU/64 GiB/bizon,
   22081_1 FAILED 1:0 on 32 CPUs/bizon and 22155 FAILED 1:0 on two CPUs/bizon.
   Compare the original failure identity with the retained preflight. Preserve
   its successful preflight/runtime equality and reviewed recovery executor
   commit/source identity. Check the original executor's relevant tree is clean.
2. Reuse the exact original plain-data completed-job and parent contracts from
   their hash-bound historical source. Execute only those two function definitions
   with the original interpreter command string explicitly supplied; do not
   import its scientific dependencies or mutate process-global `sys.executable`.
   Require three completed original phases in order, correct commands/logs,
   reused/recovered stage identities and unavailable original statistics.
3. Recheck all original/new prerequisite, source, numeric/input/graph, package/base,
   import/log/output and Git-bound file identities. Both new receipts must remain
   unadmitted observations. Require the successful single unchanged-runner command,
   original working directory/source, debug allocator/default GC/thread controls,
   helper preparation lineage and its actual child completion JSON.
4. Independently stream every saved graph array without importing NumPy or the
   native optimizer. Accept only bounded one-dimensional non-object NPY v1/v2
   arrays with little-endian int32 endpoints/float64 weights. Validate lengths,
   endpoint ranges, finite values and exact EOF. Recompute all ordered canonical
   int64 endpoint, raw int32 constructor-pair and float64 weight hashes. This is
   a narrow typed-array reader, not a general NPY parser; its format follows the
   [NumPy specification](https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html).
5. Require 984,137 vertices, 25,501,180 edges and the original three fingerprints.
   Validate exactly one completed native optimizer call, unchanged before/saved/
   after graph, CPM 0.12, seed 4, isolates, one-CPU context and imported scientific
   module/native/interpreter identities. Check the constructor's int32 input shape
   and original exhaustive parity observations, not just dimensions or pre-call
   success. Tie the resulting optimizer partition to the recovered strict seed.
6. Require both original refinement reports to match the completed new-runtime
   metadata apart from destinations. Check exact output bytes, zero selected
   directed hits, the same 984,137-gene/78-species numeric checkpoint, and full
   unique membership in all four reused/recovered partitions and the new repeated
   output. Counts are fixed to their original evidence, not new scientific endpoints.
   Recheck all bound records after validation before setting `seed_admitted=true`.

## Output And Limits

Successful status: `cpm_helper_runtime_recovered_seed_admitted_unscored`.
Record separate original optimizer and new independent-refinement identities,
all four stage origins and output hashes, complete membership checks, recomputed
graph/constructor hashes, current accounting and missing profile/timing statistics.
Preserve `accuracy_evaluated=false`, `downstream_admitted=false` and
`publication_ready=false`. Report zero native attempts by this read-only audit.

Run on one permitted Threadripper CPU with 8-GiB address space, 900 CPU seconds,
20-minute outer wall bound and disabled core dumps. Do not launch timing or
change services/unrelated processes. Metadata/file/graph checking is not an
end-to-end inference or controlled comparative resource measurement.
No DGX access, new installation or completed native diagnostic rerun.

After the audit, independently read back its report, recomputed graph hashes,
output identities and flags, then commit/push the result. A future candidate
handoff must explicitly select and bind this new admission without weakening
the historical failed-job gates. Candidate preparation/admission, phylogeny,
pair conversion and scoring remain separate. Only fully admitted high-CPM
accuracy may join the original seven-arm parameter/18-endpoint uncertainty
family; preserve seeds, multiplicity, failures and missing-result history.

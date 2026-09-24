# Recovered Search Admission Wrapper

`benchmark_tools/admit_blast_recovery_search.py` implements a separate
recovery contract; it does not modify the original search status or claim
that interrupted job 21713 completed successfully.

Before creating any validation output it requires completed merge job 22150,
exit 0:0, bizon, two CPUs and 64 GiB; the exact detached merge executor
revision `a449ff580aca58e8d1fdc653e7615e007093771a`; clean executor source;
matching merge source identities; terminal candidate-only merge status;
the expected candidate/log paths; and unchanged recorded input/output bytes.

The wrapper then reruns the frozen merge prerequisite checks without copying
the candidate: source/binary review, prefix recheck, all native/admission jobs,
complete replay panel and original partition. It compares fresh diagnostic
selection and coverage with the executed merge. The existing database auditor
extracts every formatted sequence and requires exact parity for all 984,137
corrected inputs. The integrated content auditor checks the complete recovered
table, query-block inventories and no-hit/failure dispositions. Final counts,
diagnostic signatures and before/after identities must agree.

Only successful completion sets `search_admitted=true`; accuracy, publication
readiness and downstream-execution authorization remain false. The new report
retains the candidate path rather than replacing the original partial search.
No BPO job is automatically released. Recovery limitations, logged failures,
historical durability uncertainty and uncontrolled timing remain explicit.

Ninety-nine focused recovery tests pass. New tests cover execution allocation,
source/status/output binding, malformed or duplicate accounting, existing-output
preservation and pre-output rejection of a pending merge. A real prerequisite-only
probe against current Slurm accounting rejected pending job 22150, with no output
directory created. Full wrapper orchestration/failure-injection tests and frozen
scheduled integration remain required; no production admission was attempted.

Proposed command after those integration checks, from a frozen executor:

```bash
python -m benchmark_tools.admit_blast_recovery_search \
  --root /absolute/path/to/orthohmm \
  --output /fresh/recovered_search_admission
```

## Scheduled Integration Follow-Up

The orchestration tests subsequently passed: successful sequencing plus seven
injected failures at prerequisite validation, database parity, table audit,
row reconciliation, diagnostic selection, repeated plan verification and
post-validation byte mutation. Every injected failure records a failed report
with search/downstream authorization false. These use controlled fixtures,
not production database extraction or native scheduler execution.

All 107 focused tests also pass from the detached executor at
`a21c65f2b449d828e3fefc86148cbf4fd87b6ade`. Validator **22151** is queued
with `afterok:22150`, two CPUs, 64 GiB, 24 hours and no requeue.
[Submission receipt](qfo_blast_recovery_search_admission_submission_22151.json)
records its frozen script, executor and verified pending scheduler allocation.
No production admission has executed yet, and the held BPO/inference chain
remains unchanged.

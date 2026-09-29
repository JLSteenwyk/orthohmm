# Native Completion Boundary

The collector previously waited for the launched parent and then took final
resource observations, without requiring an anchor-only native user subtree.
A parent can exit successfully while a child continues, including a child
that starts a separate process session. Parent exit alone therefore cannot
establish completion of the inference workload.

Collector v5 checks the already-retained initial and post-command cgroup
thread inventories: both must contain exactly the parked worker's thread,
with the same start-time identity, unchanged scope, valid affinity evidence,
no read errors and correctly ordered observation times. It writes
`native_completion.json` before the final memory sample and rejects incomplete
evidence. Failure cleanup releases only the owned step wrapper; the failed
attempt is retained, never silently retried or admitted as complete timing.
Replay recomputes this check from the raw points and matches both the separate
receipt and embedded record. Historical v1-v4 records remain unchanged and
are not retrospectively given this evidence.

The [kernel cgroup documentation](https://docs.kernel.org/admin-guide/cgroup-v2.html)
defines thread membership and its non-atomic read behavior. These boundary
observations do not prove continuous containment, account for work migrated
outside the subtree, or establish complete-job teardown accounting. The
existing final-job accounting gap is not removed by this change.

## Owned Diagnostic

One 1-CPU, 128-MiB, two-minute-limit diagnostic ran on `bizon`, job **22372**.
Allocation and batch step both COMPLETED 0:0, elapsed one second; subsequent
queue query was empty. The [raw result](threadripper_native_completion_22372.json)
retains commands, source hashes, anchor/child identities and boundary reads.
A `/bin/true` clean control passed. A Python parent then launched a 20-second
sleep in a new session and exited zero. The live child appeared in the native
subtree inventory, and the completion check rejected it. The diagnostic
terminated only that identified owned child through its pidfd.

The [executed diagnostic source](probe_native_completion_22372.py) preserves
the exact script bytes pinned in the result. The working diagnostic was then
hardened to avoid signalling if owned-scope validation itself fails; it was
not rerun. No failure was encountered on that path in job 22372. The diagnostic
did not run the full collector, native OrthoHMM or OrthoFinder, and is not a
production timing repeat or measurement-overhead test.

146 focused tests pass across collection, replay, thread observations,
disk-backed history, executor, resource wrapper and reporting. Cases include
live descendants, changed/missing identity or scope, nonfinite/reversed times,
missing/altered/symlinked receipts, failed-attempt retention and cleanup.
Full-scale/native integration, final accounting, the reviewed environmental
policy and a verified quiet window still precede the 27 production runs.

# High-Sensitivity Four-Proteome Third Repeat Submitted

After the [canonical index-18 review](THREADRIPPER_SHARED_ATTEMPT_22414.md)
passes all four categories and its reporting milestone is committed/pushed at
`cf6323ab`, the unchanged new continuation helper independently checks the
complete nineteen-attempt history and releases exactly index **19** as job
**22415**. This is high-sensitivity OrthoHMM, four proteomes, repeat **2**,
73,266 frozen input proteins. Neither earlier excluded attempt is retried.

The [public submission](threadripper_shared_submission_22415.json) exactly
matches canonical launch metadata: 1,179 bytes, SHA-256
`6a16b73d9a992146e8fee7d2797645a55c30db26c16d4322d96a332be7b81c98`.
Actual held-job validation, request selection, scheduler-comment digest binding
and release complete successfully. Request is 6,881 bytes, SHA-256
`8306dc3e6ccae64d7983787b5d70f19f67bc0b40f6ece58d7124adf63309d3e3`.
The same scientific configuration, repaired runtime/source/resource/readiness
binding, planned order, 64-SMT-slot allocation for 32 physical native CPUs,
128 GiB and 26-hour scheduler limit are retained. No exclusive allocation,
requeue, displacement or change to unrelated work.

## Actual Native Handoff

The [live snapshot](threadripper_shared_live_22415.json) is 12,006 bytes,
SHA-256 `5f81a75a0679c74cb69452fe941005ba16c401bce40477a122bb4cac6b60746d`.
It contains a fresh RUNNING controller observation at 3:33, verifies the real
request digest and live native-worker job/step membership, and records affinity
0..31. Native log confirms built-in profile-HMM/k-mer search, high sensitivity,
BLOSUM62, E value 0.0001, Leiden CPM 0.1, 32 CPUs and `--stop infer` without
phylogenetic inference. All-to-all comparisons are running, not completed.

Actual environmental preflight passes with available memory
329,388,662,784 / 329,634,897,920 bytes, above the frozen 137,438,953,472-byte
floor. Foreign demand is 41.1505908056 CPU-core equivalents. Release binds the
preflight and `go` is true. The immutable initial process observation remains
pinned while periodic data append; first observed start gap is 30.0004490131s,
within the unchanged 35-second bound. This is early repaired handoff evidence,
not whole-run cadence/resource/output review, isolation or a slowdown estimate.

The current resource table stays at nineteen reviewed identities; this live
index 19 is not a reviewed timing. Obtain fresh scheduler state on continuation.
Wait for parent/batch/native terminal outcomes, preserve terminal controller
evidence before purge if needed, then use the new independent reviewer on index
19. Do not resubmit, restart after an observation timeout or advance index 20
before canonical review and full-prefix validation pass. All owned tool sessions
are terminal; Slurm continues independently. Full publication work is incomplete.

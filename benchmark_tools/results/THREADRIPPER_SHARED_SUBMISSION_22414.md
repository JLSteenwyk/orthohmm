# Shared Panel Resumed: OrthoFinder Four-Proteome Repeat

After the [tested pre-native repair](THREADRIPPER_PRENATIVE_REPAIR_20261004.md)
and [refreshed runtime/reporting integration](THREADRIPPER_PRENATIVE_CONTINUATION_20261004.md),
actual preparation and independent binding/policy readback pass. The public
[preparation](threadripper_shared_prenative_preparation_20261004.json) exactly
matches canonical bytes. All 887 execution sources are pinned at `af32242d`;
the native startup/input/scientific configuration is unchanged. Original
index-0 and index-17 failures remain excluded; neither is retried.

The [submission receipt](threadripper_shared_submission_22414.json) exactly
matches canonical launch metadata: 1,176 bytes, SHA-256
`0878bbc9abba9937952e1638772f5e04ddc2678cc4737729d12ae40a9180734e`.
Job **22414** is index **18**, full OrthoFinder **3.1.5**, **four proteomes**,
repeat **2**, 73,266 frozen input proteins. The held job is independently
selected against the complete eighteen-attempt history, then bound to request
digest `97345039e3cf94d45d89d48a88ad18cea3cb1e38de3a54a535be320e988840f9`
in its scheduler comment before release. It has the unchanged 64-SMT-slot
allocation for 32 physical native CPUs, 128 GiB and 26-hour scheduler limit,
without exclusive scheduling, requeue or changes to unrelated work.

## Actual Native Handoff

The [live snapshot](threadripper_shared_live_22414.json) is 12,598 bytes,
SHA-256 `ecd90b157e0e0d5f532b6e2b40ee4a6c67c551f778246a9f44b18d00b771eb15`.
Its fresh scheduler query observes RUNNING at 3:38 and validates the actual
request digest. The selected native worker is live in the proper job/step
cgroup with affinity 0..31. Native log confirms OrthoFinder v3.1.5, 32 threads,
default MSA tree inference, successful MCL/FAMSA/FastTree program checks and
DIAMOND all-versus-all search. This is full inference, not the sequence-only
pre-phylogenetic checkpoint.

Actual preflight passes and records available memory of 329,769,832,448 /
329,574,883,328 bytes, above the prospective 137,438,953,472-byte floor.
Foreign demand is 41.2694294498 CPU-core equivalents. Environmental release
binds that preflight and `go.json` records `go: true`, not an abort. The
immutable initial process observation is pinned by the snapshot while the
process stream continues growing. Its first observed start gap is
30.0004067310 seconds, below the unchanged 35-second bound. This establishes
the real repaired handoff, not whole-run cadence or final resource validity.

The snapshot is dated live evidence, not a permanent scheduler-state claim.
Obtain a fresh query on continuation. No terminal resource/output review or
final timing is available yet. Wait for parent/batch and native completion,
then use `review_shared_prenative_panel_20261004.py --index 18` with the frozen
private controller. Capture terminal controller evidence before scheduler
purge if necessary, and independently corroborate retained evidence with
fresh accounting. Do not resubmit, retry index 17, invoke the old launcher or
advance index 19 before the canonical review passes.

The current table remains the dated 18-reviewed-attempt checkpoint; index 18
is not yet reviewed. Shared-host distortion is unknown and potentially
method-dependent. No timing adjustment, causal speedup, isolated-efficiency
ranking or publication-completion claim follows. Remaining panel, final
manuscript reconciliation and versioned/archive release remain active.

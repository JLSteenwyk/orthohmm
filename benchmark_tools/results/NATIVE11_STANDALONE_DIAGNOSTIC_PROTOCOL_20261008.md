# Native11 Standalone Semantic Diagnostic

## Failure And Question

Native23985 completed0:0 in17:28:07. Original reviewer23986 failed0:11 in
00:45:50 without review.json or failure.json. Its retained faulthandler output
locates garbage collection during the materialized/final-clustering partition
comparison at validate_native_factorial_outputs.py:331. Runtime, scheduler,
resource replay, resources and environment intermediates are retained, not
accepted as a completed review. Location and loaded extension names do not
establish cause. No OOM, malformed partition or runtime defect is proven.

The question is whether the unchanged allocated output validator completes in
a separate process against these same retained outputs. Historical22734
standalone semantic success for a different cell is not evidence for native11.

## Prospective Execution

[Batch](native11_standalone_diagnostic_20261008_v1.sh) runs only the unchanged
allocated semantic validator production CLI. Use the retained scientific
Python3.10 environment, faulthandler, ordinary garbage collection, sanitized
Python/loader environment and one thread per scientific backend. No interpreter
upgrade, dependency installation, GC suppression, monkeypatch, skipped check,
new native inference or full resource/environment replay is authorized here.

One new held Slurm job:2CPUs/32GiB/6hours, bizon/gpu, no array or requeue.
Prospective sources/tests must be committed and pushed before submission.
Record actual submission and inspect owner, command, working directory,
comment, resources, absent dependencies and fresh output namespace. Commit
and push the held receipt before one release; recheck safe capacity and bound
sources/request first. Ordinary contention is accepted and recorded, not an
exclusion or reason to require the DGX. No unrelated job/service changes.

Pinned request:
`benchmarks/work/native_factorial_launch_20261004/request_11_allocated_v1.json`
SHA256:`7bf63b80bd5932b9edbd1b2c5ff3fb77f5557e4f6c64e077045d6a50c8d366a1`.
Pinned allocated validator SHA256:
`146f1333715b5ef52ee17de855fad584d1efc91c7a1154edcd269829b5a871be`.
Pinned unchanged semantic kernel SHA256:
`3357637503f35238f654edea5c4c12bd293f8d20908c141421a2d40278af7c8d`.
Pinned original reviewer SHA256:
`c8707fdb1855a0d73fba537c001b4a355faff094a38a177a05e179148924c929`.

The production validator repeats its request/amendment/plan/helper/input and
successful-native-terminal gates, placement and frozen native settings,
checkpoint contents, universe, partitions and candidate-sidecar checks.
Write only a fresh diagnostic namespace:
`benchmarks/work/native11_standalone_diagnostic_20261008_v1/`.
Retain GNU time, stdout/stderr, scheduler outcome and actual output if produced.
Retain every original23986 file; do not overwrite or relabel that failed attempt.

## Interpretation And Recovery Boundary

If this diagnostic fails, retain its actual stack or exception and use that
evidence for the next focused diagnosis. No automatic retry. If it succeeds,
conclude only that standalone validation passed in this invocation. Success
does not explain the original crash, prove a fix, validate a full review,
admit accuracy/timing, or authorize native12. Any subsequent full-review
recovery requires its own prospective justification, unchanged original gates
and fresh output namespace; this protocol does not submit such a review.

All diagnostic timings describe shared-host postprocessing, not native
inference cost or isolated efficiency. Competing CPU, memory-bandwidth and I/O
work may have an unknown, tool-dependent impact. Full publication requirements
remain active and completion unproven.

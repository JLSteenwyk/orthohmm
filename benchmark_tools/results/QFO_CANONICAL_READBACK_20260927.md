# Canonical QfO Readback Integration

Added `readback_qfo_canonical_phylogeny.py` to connect canonical scheduler and
artifact admission, four independent scientific readers, and the prespecified
three-way native prediction comparison. It takes the actual future job ID,
plan SHA and frozen submission SHA; no canonical job ID is invented here.

Admission verifies the completed 32-CPU attempt, source/plan/environment/rules,
reuse-admission identity, exact outputs and complete input inventory. Copied
raw-tree artifacts must come from the admitted fresh run with unchanged tools.
The pre-execution copy receipt is verified against original bytes; final
checkpoint bytes are validated separately because reconciliation rewrites
checkpoints. Species-tree and reconciliation-output reuse are rejected.

The four readers validate canonical output. Full historical/fresh/canonical
contrasts then report root partitions, affected root and pair families,
native pair/confidence differences, species-tree topology/bytes, summaries
and membership audits. Runtime artifacts are rechecked after the readback.
This does not score benchmark references or replace historical metrics.

Validation: 53 focused tests pass. A
[real fixture receipt check](qfo_canonical_cache_receipt_fixture_20260927.json)
verifies the four actual copied files against the prior fresh fixture outputs.
This is not full canonical runtime admission; the scheduler-integrated path
awaits the real canonical job and complete outputs.

Native 22329 remains RUNNING; readback 22330 is pending on its dependency.
No canonical plan/job exists yet. The readback sources for those jobs were
left unchanged. Next use successful full fresh readback to prepare and freeze
the canonical plan, submit it once, then freeze this reader's dependencies
and exact submission binding for its scheduled audit. DGX remains deferred.

# Corrected Native OrthoMCL Inference Queued

Submitted **21750**, `afterok:21749`, requesting 180 CPUs, 900 GiB RAM,
seven days on bizon, with no requeue. Authoritative `scontrol` state is
PENDING/Dependency with zero restarts. No corrected native inference has
executed yet. The resource request matches the prepared180-thread setting;
the native inter-species phase uses64pair workers, not180parallel MCL jobs.

Frozen runner: `benchmarks/work/publication_qfo_corrected_orthomcl_v1`,
revision `22f4eec4314304a0cda71e02595adfb20827c34b`.
Batch: `benchmark_tools/results/qfo_corrected_orthomcl_infer_batch_20260918.sh`.
BPO admission job argument:21749. Its frozen validator remains
`914c3fe6f8618c3d90452e7a5b6d7f03f939ccd2`.

## Execution Contract

`run_qfo_corrected_orthomcl.py` requires completed independent BPO admission,
its fixed executor/source, scheduler provenance and exact admitted BPO,
offset/query indexes and original GG records. It checks the corrected input
scope, frozen configured native sources, clean tool directory, dedicated
Python and native Perl/system runtimes, and space for staged inputs plus a
10-GiB reserve. This reserve is not a prediction of total graph-output size.

It creates separate verified input copies in a fresh inference directory,
checks all copied indexes, and invokes guarded native mode4. It does not
rerun BLAST or switch scientific defaults. The prepared source manifest is
SHA-256 `9ccd23429301d754f80fc72ae4ba4b28ef5036116e0670c5aebdb85a41b00ae6`.
Existing inference outputs/caches are not silently resumed.

Native execution is timed separately with `/usr/bin/time -v`. The wrapper
checks the exit code, requires one final group file, rejects malformed,
unknown or repeated memberships and duplicate labels, and reports ungrouped
inputs rather than silently adding them. It requires all3,003nonempty
species-pair cache files for78species and retains their hashes. Input/index
hashes and modification times, source records and runtimes are rechecked
after execution. Failures remain recorded.

Expected execution report:
`benchmarks/results/qfo_corrected_orthomcl_v1/inference_execution/status.json`.
Success remains `corrected_orthomcl_native_exited_zero_pending_admission`.
Native groups under the configured tool directory are not yet admitted for
scoring. Source-query failure diagnostics remain attached to the report.

## Verification and Boundaries

The full unit suite passed: **3,859 tests in85.70s**, with installed legacy
runtime tests enabled. New runner tests cover admission/source/hash/scope
binding, stale directories, allocation, worker count, insufficient space,
group memberships, cache inventory and retained execution failures. Actual
native components remain supported by the separately retained staged-input
and serial/parallel fixtures; orchestration tests use mocks.

Batch syntax and clean frozen-worktree checks passed. The frozen runner's
check-only invocation against still-pending21749 correctly rejected
nonterminal evidence and left no inference directory. This is a fail-closed
gate test, not a successful full-data preflight.

Next: independently admit terminal native outputs, verify final MCL group
semantics and conversion, score the corrected reference, and reassess the
impact of logged BLAST failures. Shared-host stage accounting is not matched
end-to-end timing; small-fixture agreement does not prove all64-worker
schedules or biological superiority. Publication readiness is not claimed.

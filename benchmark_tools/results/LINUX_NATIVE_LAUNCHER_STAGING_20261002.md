# Native Constructor Test Launcher Staging

The first actual [Linux diagnostic CI](LINUX_NATIVE_DIAGNOSTIC_CI_20261002.md)
job, source-b7ed/run 36980170796/job 110752659858, fails with 75 passes,
one failure and zero skips in 13.50s. The frozen-import high-CPM child exits
1; the original test hides captured stderr in a CalledProcessError. Dependency
installation succeeds and the failed JUnit artifact is retained. Do not
interpret this as a successful 76-case runner or restart it.

The test points at an ignored local launcher, not a committed checkout file.
A bounded local missing-launcher reproduction exits 1 with `Wrong frozen
worker import` and creates no result. This establishes that failure mechanism
locally, not the original child's exact hidden traceback. Seven constructor
and twelve repeat-worker cases passed remotely; preserve the one failed case.

Stage the same five committed package source files already used by the other
native worker tests into a temporary launcher. Run the unchanged diagnostic
with that explicit root, retaining actual igraph/Leiden imports, affinity,
loaded-library records, graph fingerprints and all no-optimizer/accuracy
checks. Assert the frozen worker module comes from the staged file. Do not
create the old workstation directory, reuse its runtime, fabricate a worker
stub or change production admission. Expose captured stderr through the
return-code assertion for future failures.

**76 local cases pass in 20.05s**, zero failures/errors/skips, including all
twenty native worker variants. This follows the first local 317-case pass;
it is not another independent scientific experiment or full regression.
The existing no-skips CI gate and workflow remain unchanged. New-runner
confirmation of this staging fix remains pending.

[Source/report/log identities](linux_native_launcher_staging_20261002.json)
retain the remote failure, local negative and final test receipt. Only the
test and documentation change relative to b7ed. Production diagnostic/core,
scientific settings/scores, historical archives and receipts stay unchanged.
Shared installed dependencies are not a hermetic or full frozen-runtime
restoration; no full-scale stability or complete release claim follows.

The preceding macOS 241-case confirmation remains valid for its source and
scope. Other workstation/raw/platform failures, other-QfO uncertainty,
dependency/data rights, controlled resources, full release and public
deposition remain open. Timing stays deferred; no contention poll, scheduling
question, DGX, unrelated process/service action, biological inference, raw
scoring, bootstrap, archive regeneration or source search occurs.

# Host Reader Platform Coverage

The live host-reader test requires Linux `/proc` and unified cgroups. Its
failure on macOS is a test-placement error, not a portable counter interface.
Mark only this live integration test Linux-specific and add the complete
host-reader module to the existing Ubuntu Python 3.12 diagnostic lane.
Keep all parser, summary and synthetic-reader tests on both configured platforms.
Production reader code and historical resource/scientific evidence are unchanged.

The Linux JUnit gate requires exactly 35 host-reader cases, including exactly
one live probe/no-overwrite case. It also retains the existing eight native
constructor and 12 native repeat-worker requirements. No selected case may
be skipped, failed or errored. The full lane should execute 111 cases; a missing
Linux capability fails rather than quietly skipping the live reader. This
does not remove the four external raw-data exporter failures or provision data.

## Verification

Add 17 synthetic cases: snapshot/raw/optional/read-bracket preservation,
missing or inaccessible optional counters, missing required files, changed
or malformed cgroup membership, invalid CPU evidence, invalid intervals,
synthetic run serialization/no-overwrite and source mutation. Missing optional
values remain absent with explicit errors, never zero-filled. Required failures,
changed scope and source mutation prevent a successful receipt.

**110 local cases pass in 21.36s**, zero failures/errors/skips, with one
deliberately deselected live host test. This includes the five existing native
diagnostic modules and 34 synthetic/parser host-reader cases. The earlier
34-case pass in 0.20s overlaps this panel. Do not count the deselected test
as a pass or infer 111 local executions.

An isolated copied Python 3.10.13 child passes all 34 host synthetic/parser
cases in 0.08s, again with the live case deselected. Seven files/loaded project
origins are staged. A Python open-event guard blocks three canaries (original
checkout, `/proc/stat` and `/sys/devices/system/cpu/online`), observes zero
later forbidden opens, and prohibits child subprocesses. Temporary staging
is removed. The actual decorator expression evaluates to no skip on Linux and
skip on macOS/Windows; this expression check does not emulate those kernels.
The guard is Python-event evidence, not OS containment or native host validation.

```bash
python -m pytest -q tests/unit/test_probe_host_counters.py -k 'not live_read_only'
```

[Source, workflow and receipt identities](host_reader_platform_coverage_20261002.json)
pin this evidence. The live Linux execution is a CI requirement, not yet a
verified result for this patch. Do not run the live test on the Threadripper
while timing/host contention polling is deferred.

## Remaining Scope

Observe preceding source-4ccd72a6 run 36994103659 without restarting it:
docs/Linux diagnostics succeed, wheel initially live, full initially queued
and all four macOS fast jobs live. Later the wheel succeeds and all five
macOS jobs remain live. That source contains the preceding raw relocation
patch, not this platform amendment. No new-patch remote confirmation, sibling
test counts or full CI success follows from these statuses.

At 10:21:56 UTC, macOS Python 3.11 is terminal/failure; the four other test
jobs remain live and docs/Linux/wheel succeed. Inspect its actual job
110796750061 log once at 10:22:27 UTC, confirming checkout SHA: 14,007 passes,
five failures, 118 skips, zero errors and 30 warnings in 390.43s. All 43
preceding raw-relocation tests pass, as do all 44 inventory/fixture cases;
do not infer all 187 local raw-panel cases remotely. Four default raw exports
and the Linux host-reader case fail. The new platform amendment is not in
that source, and the synthetic binding tests are not native raw restoration.

This live read-only availability check on an ephemeral Linux CI runner cannot
establish the Threadripper's quiet window, observer overhead, controlled timing
or restored scientific runtime. All production admission flags remain false;
no timing, service/scheduler change, DGX access or unrelated job action occurs.
Raw data provisioning/rights, other-QfO uncertainty, controlled resource evidence,
complete release and deposition remain open. The publication goal is incomplete.

# Additional Launcher Test Environment Isolation

Two more synthetic test modules call standalone launchers that intentionally
change process environment. An explicitly requested cleanup fixture now restores
that state for their six parameterized launch cases. Production launchers,
historical plans, native commands, scientific settings and scores are unchanged.
The fixture is not autouse and does not mask failures in unrelated tests.

## New Evidence

The existing source-c98 macOS Python 3.11 fast log contains **13,718 passes,
47 failures, zero errors, 111 skips and 30 warnings in 361.60 seconds**. All
nine replay cases pass, confirming the prior temporary-script relocation.
The earlier 94-case scope gives 81 passes, 12 failures and one skip: five
malformed-invocation and seven Pandoc cases still fail; the separate baseline
wheel compilation case skips. Its skip cause is not established. Setup logs
GCC 16.2.0, GNU Time 1.10 and Pandoc 3.11, but that alone does not prove later
availability or complete-suite success. Only this fast-job log is inspected.

An isolated local child observes the six pressure/frontier-overhead launch
cases before the new fix. All six pass, yet each leaves environment deltas.
The first changes nine keys, including PATH, Python configuration and thread
limits; the remaining five change the pycache prefix. This is passing-but-leaking
test behavior, not a new biological result. All measurement and working-directory
calls remain mocks; there is no remote access or host workload observation.

Pressure cases occur before the failing Pandoc/Bash-dependent cases in the
actual CI order. Their launcher replaces PATH with historical tool paths plus
`/usr/bin:/bin`, omitting Homebrew directories. That is a plausible mechanism
for the remote symptoms, but the old log does not contain an environment
readback after those tests. Causation and compiler-skip diagnosis remain
inferences requiring fresh remote confirmation; sibling causes are not assumed.

At 04:38:19 UTC all five source-c98 test jobs have failed and wheel/docs succeed.
Source-deed still has five live test jobs and successful wheel/docs. No existing
run is restarted, no sibling test log is downloaded and no live job is treated
as terminal because of an observation timeout.

## Correction And Validation

`isolated_launcher_environment` snapshots the complete original environment and
restores it in `finally`. Only the two mutating test functions request it.
Four new guards test normal/exceptional exits and both relative teardown orders
with an inner monkeypatch, checking full-map equality for added/deleted/changed
keys. The prior frontier-native fixture's separately tested cleanup stays intact.
Existing recipe, authorization, argv, failed-command and retained-output checks
are not relaxed, skipped or replaced with new scientific admission.

The expanded ordered local panel gives **183 passes in 8.43 seconds**, zero
failures/errors/skips and no after-test environment deltas. Module counts are
frontier-native 17, pressure 47, frontier-overhead 38, cleanup guards four,
replay nine, malformed-invocation eight, manuscript 13, bibliography 19,
relocation helper eight and CPU-wheel 20. Real local Pandoc rendering and
baseline compilation execute after all three launcher modules. The observer
does not restore state and excludes only pytest's bookkeeping key.

These observations include the earlier 94-case scope; panels are not additive.
They are local, not actual macOS confirmation of the new correction or a full
regression. The [receipt](ci_additional_launcher_environment_20261002.json)
pins both JUnit reports, test sources, unchanged production/plans and the actual
remote log. The pre-fix passing-but-leaking evidence is retained.

## Remaining Work

Inspect the next automatic CI log for this expanded scope before claiming remote
closure. Remaining workstation snapshot-path/raw-data, Linux-specific capability,
and complete executable-release failures are still open. The 55-file direct-review
archive remains separately verified; it is not transitive study evidence or
release clearance. Other-QfO uncertainty, original TreeFam sources, dependency/
data rights and controlled resources remain incomplete. Timing stays deferred
without contention polls, quiet-window questions, DGX access or actions on
unrelated jobs/services. The full publication goal remains active.

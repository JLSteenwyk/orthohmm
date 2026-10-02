# Launcher Test Environment Isolation

The in-process historical launcher tests passed while leaking process environment
changes into later tests. Restore the complete original environment in the
fixture's `finally` block. This test-only correction preserves production
launchers, archived scripts, scientific settings, scores and admission rules.

## Observed Behavior

A bounded child-process probe of the three launcher cases passes before the fix
but observes environment changes after each test. The first changes nine keys,
including PATH, Python configuration and thread limits. The second changes the
pycache prefix; the third changes PATH and that prefix. The observer records key
names and PATH hashes, not environment values, and does not perform restoration.
Only pytest's own `PYTEST_CURRENT_TEST` bookkeeping is excluded.

The synthetic launcher recipes replace PATH with historical tool directories
and `/usr/bin:/bin`. Losing the macOS Homebrew paths is a plausible explanation
for later BSD Bash, missing Pandoc and compiler symptoms. This is an inference,
not a remotely demonstrated causal result: the old CI log does not include an
after-test environment audit. No remote compiler-skip cause is established.

At 04:20:06 UTC on 2 October, source-1ce CI has five failed test jobs and
successful CPU-wheel/docs jobs. Its actual Python 3.10 fast log contains
13,715 passes, 48 failures, zero errors, 111 skips and 30 warnings in 458.40s.
GNU Time 1.10, Pandoc 3.11 and GCC 16.2.0 are logged during setup. The previously
changed 158-case utility/reference scope now gives 152 passes and six explicitly
retained-data skips, without failures/errors. That confirms the prior fixtures
on actual macOS, not this new fix, full-suite success or biological inference.
Only that fast-job log is inspected; sibling causes are not inferred.

## Focused Corrections

The fixture restores added, removed and changed environment keys on normal and
exceptional exit. Two new cases check exact whole-environment equality before
and after the inner monkeypatch teardown. Existing synthetic host, measurement
and working-directory mocks remain; no DGX access or host workload poll occurs.

The corrected-replay batch test now uses the existing explicit relocation
helper with a temporary root containing spaces, a dollar sign and an apostrophe.
Only a temporary copy's root/interpreter bindings change. The original archived
script stays byte-identical, and successful argv must contain the exact root.
Commit/hash/dirty-executor failure guards and the no-inference stub stay intact.

An existing helper test assumed the current interpreter could not equal the
historical interpreter. That is false on this workstation. It now compares the
entire expected rewritten script, retaining the two-call and source-immutability
assertions; the helper itself is unchanged.

## Executed Checks

The ordered local panel runs launcher, replay, malformed-invocation, manuscript,
bibliography, relocation-helper and CPU-wheel tests. **94 pass in 7.53 seconds**,
zero failures/errors/skips and no after-test environment deltas. Module case
counts are 17, nine, eight, 13, 19, eight and 20 respectively. Actual local Pandoc
rendering and baseline compilation execute after the launcher. These are local
test observations, not a fresh macOS or complete-suite result.

The initial ordered panel gives 93 passes and one failed interpreter-assumption
assertion, already with no environment deltas. Its report and the pre-fix
three-case probe remain retained. Panels overlap and must not be added.
The [machine-readable receipt](ci_launcher_environment_20261002.json) pins all
three JUnit reports, changed test sources, unchanged production/helper/archive
sources and the inspected remote log. JUnit suite durations and console
durations differ slightly; the receipt preserves both rather than rewriting XML.

## Remaining Work

Observe the next automatic CI without restarting existing runs. Full regression,
versioned executable release, source/rights, uncertainty and comparable resource
evidence remain open. Controlled timing stays deferred while other publication
work proceeds, without quiet-window questions, host polling, shared upgrades or
actions on unrelated jobs/services. No publication-readiness claim follows.

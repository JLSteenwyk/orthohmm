# Portable QfO Batch Handoff Tests

Five QfO handoff-test modules launched dummy executors through a historical
workstation Python path. Relocate only that interpreter and `ROOT` in their
temporary copies; preserve the archived batch scripts and all handoff guards.
This improves executable test portability, not production workflow portability
or the scientific evidence for OrthoHMM.

## Evidence And Correction

The retained source-df22 macOS Python 3.11 log contains 38 failures across these
five modules, each reporting missing `/home/bizon/anaconda3/bin/python`. Reuse
that already downloaded log, rather than restarting or downloading old jobs.
At the 02:07:15 UTC snapshot, the subsequent source-57a8391f CI run 36953675242
has successful CPU-wheel/docs jobs and five live test jobs. That source includes
the prior boot-fixture corrections, not this new handoff correction; outcomes
and causes of live jobs are not inferred.

An explicitly requested test fixture writes a fresh copy, checks the expected
root/interpreter markers, quotes the replacement paths and refuses overwrites.
The two Python invocations in reconciliation admission must both be relocated.
Its fake `sacct` command also uses the test interpreter through a quoted shell
launcher, not a workstation shebang or native scheduler call. Executors and
roots now contain spaces, dollar signs and apostrophes in all handoff cases.
Expected argv, input hashes, thread limits, commit/dirty-source checks, task
identity, index/parity rules and failed/running/duplicate accounting rejection
remain unchanged. Only dummy executors and synthetic accounting are exercised.

## Executed Validation And Boundaries

**329 focused cases pass in 14.95 seconds**, zero failures/errors/skips.
This includes all 238 cases in the five handoff modules, eight new fixture
checks and 83 prior clock/accounting cases. The earlier overlapping 248-case
panel is superseded, not additive. New checks cover exact permitted edits,
shell execution with special root/interpreter paths, missing/duplicate markers,
overwrite refusal and explicit two-call relocation. The shared boot fixture
and pytest configuration retain identical ASTs; all eight historical scripts
remain byte-equal to the prior commit. The
[receipt](ci_qfo_batch_fixtures_20261001.json) pins sources, scripts and JUnit.

Commit/push this tested milestone and observe its own automatic CI handle,
without restarting older runs. New remote handoff execution remains pending.
Historical production batches still require their recorded host/runtime; this
test-only fixture is not an authorized production relocation tool. Final release
workflow portability, broader CI, raw/native runtime and rights closure,
uncertainty/source gaps, controlled timing and archival deposition remain open.
No benchmark/scorer/engine/default changes or new native runs occur here.
The Threadripper timing panel remains deferred without host polls, quiet-window
questions, DGX access or unrelated workload/service modifications.

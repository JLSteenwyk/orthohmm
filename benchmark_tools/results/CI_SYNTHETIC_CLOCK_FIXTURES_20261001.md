# Synthetic Accounting Clock Fixtures

Continue from pushed `df22e79ac82cd67feb43343c84ee28f28e2f42f2`.
The user does not know when the Threadripper will be quiet and asks whether
that is needed now: no. Defer controlled timing and continue bounded
reproducibility work. No contention poll, quiet-window question, unrelated
job/service action, shared-package change or DGX access is needed here.

## Observed Remote Evidence

Inspect source-df22 CI run 36952199808, job 110667364838 (macOS Python 3.11).
The downloaded log has **13,450 passed, 185 failed, 38 errors, 96 skipped,
30 warnings in 415.42 seconds**. All 61 portable-rendering cases, seven
independent count-reproduction cases and 20 declared-dependency cases pass.
This remotely confirms those bounded tests, not the broader suite or native
scientific reproduction. Do not add these counts to earlier overlapping panels.

Four synthetic accounting modules still read the host's Linux boot-ID file.
Their saved summaries show 36 failures/errors in total; three are downstream
missing-output/partial-evidence assertion failures rather than direct boot-file
exceptions. The source log is retained and pinned in the
[receipt](ci_synthetic_clock_fixtures_20261001.json). At 01:58:42 UTC all five
test jobs are terminal failure; CPU-wheel and docs succeed. Only the Python
3.11 job log is inspected here; do not infer sibling causes or restart them.

## Scoped Correction And Validation

Opt four unit modules into the existing `synthetic_linux_boot_id` fixture.
They already use fake cgroup counters, temporary trees or fake service replies;
their boot identity should be synthetic too. Ordinary file reads remain real.
Assert the fixed identity in successful measurement, saved-evidence replay,
frontier comparison and fake-service capture. Add three failure cases:
frontier boot-file absence preserves an invalid snapshot without counter claims;
fake-service reboot remains unapproved; missing boot fails before service reads.

**83 focused cases pass in 2.36 seconds**, zero failures/errors/skips, including
the three new cases and existing timeout, altered-evidence, scope-change and
secret-exclusion checks. The measurement tests launch only their own tiny
subprocess fixtures. Service calls are injected, not native service queries.
The shared fixture and all four production collectors remain byte-equal to
the preceding commit; production Linux capabilities and timing gates are not
weakened. JUnit and changed/unchanged source pins are in the receipt.

Commit/push this tested milestone, then observe its actual automatic CI handle
without manual retries. New remote execution of these clock corrections is
pending. Remaining affinity, native-tool, raw-input and historical-path failures
require separate diagnosis. No quiet-host certification, benchmark admission,
new scores, defaults, controlled timing or publication-readiness claim follows.
The full goal stays active; source/rights, uncertainty, comparable resources
and the final public release remain unfinished.

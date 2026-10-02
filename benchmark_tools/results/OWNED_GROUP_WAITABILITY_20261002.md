# Owned Group Exit Waitability

## Observed Failure

The actual source-657 macOS Python 3.10 fast log has an EPERM failure in
`test_timeout_stops_owned_command`. The first error occurs in a zero-signal
group probe after TERM; `Popen.returncode` remains `None`. Exception cleanup
then attempts TERM again and receives EPERM with the same pending leader.
The existing helper's single `poll()` does not resolve this state. The log
does not establish the exact kernel state or exclude genuine permission denial.

The earlier [timeout correction](CI_OWNED_GROUP_TIMEOUT_20261002.md) remains
historical evidence, including its passing scope. It is not proof that every
macOS exit race is resolved. Preserve the source-657 failure and its log pin.

## Bounded Correction

On a permission error only, poll the owned leader for at most a nominal
0.1-second window, with 0.01-second sleeps. If exit becomes waitable, reap via
`poll()` and require a fresh zero-signal probe to report the group absent.
Only then accept cleanup. A leader that stays live, a still-present group or
a denied fresh probe still raises the permission error. Do not equate leader
exit with group exit, scan other processes, target individual descendants or
silently ignore errors. Normal TERM/grace/KILL handling is unchanged.

This handles the tested delayed-waitability condition conservatively; it is
not proof that this was the exact CI kernel mechanism. The timeout is a polling
deadline, not an OS scheduling guarantee. Races longer than the bounded window
and persistent denied/zombie-descendant groups remain fail-closed.

Both measurement and WGD wrappers already use this shared helper and record its
source identity. Their code is unchanged. Prospective execution inventories
must pin the new helper; no historical source/plan/receipt is repinned and no
diagnostic or scientific job is relaunched.

## Validation

**135 local cases pass in 8.86s**, zero failures/errors/skips: 15 new helper,
37 measurement, eight WGD, ten measurement-audit, 33 FAS and 32 provider cases.
The new cases cover delayed and near-bound readiness at TERM/probe/KILL, a
leader that stays alive and groups that remain present or inaccessible.
The existing native TERM-ignoring descendant test passes. Both wrapper CLI
help commands succeed. Adjacent FAS/provider checks show no test-state leakage.

Before the fix, all 15 new cases fail against the old helper: six intended
cleanup-success cases raise premature EPERM; nine denial cases retain the
correct rejection but fail the new bounded-rechecking assertions. This does
not claim nine unsafe old cleanups. The failed JUnit report is retained.
The intermediate 52-case pass and final 135-case panel overlap, not independent
replications. These tests are not resource benchmarks or broad macOS proof.

[Machine-readable source/test/report/log pins](owned_group_waitability_20261002.json)
record the evidence. At 07:00:57 UTC source-055 automatic run 36975890708 has
five live test jobs; wheel/docs succeed. No resubmission or inference from an
incomplete observation follows. The new helper is not in that source.

Later inspect only source-055 macOS Python 3.13 fast log: all previous 65
FAS/provider cases pass, including all 33 exact native-source FAS comparisons.
Overall it has 13,856 passes, 37 failures, 110 skips and 30 warnings in
477.20s. A new Python 3.13 oversized-tar fixture construction failure is
retained for follow-up; other platform/raw/retained-path failures remain.
This is prior FAS test confirmation, not the new helper's macOS validation or
a complete successful suite. Different Python-version failure totals cannot
be used as a causal regression comparison.

## Remaining Work

New macOS confirmation, other raw/workstation/platform CI failures and complete
executable release remain open. Frozen scientific source/settings/scores and
all historical results/archives are unchanged. No biological inference,
bootstrap, annotation, scoring or controlled timing is rerun. Timing remains
deferred; no host poll/question, DGX access or unrelated workload/service action
occurs. Other-QfO uncertainty, source/data rights and public deposition are not
resolved by this helper correction.

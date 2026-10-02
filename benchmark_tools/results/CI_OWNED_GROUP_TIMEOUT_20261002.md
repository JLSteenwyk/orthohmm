# Owned Process Group Timeout Cleanup

## Evidence And Change

The retained source-2d Python 3.11 macOS fast log shows `killpg(pid, 0)`
raising `PermissionError` after TERM in the measurement wrapper, followed by
another cleanup error. Its WGD timeout fails on KILL. The actual source-11d
Python 3.11 fast log independently repeats the WGD KILL failure, while all
21 measurement cases pass. Timeout behavior is not consistently failing in
every run. Both downloaded logs and local reports are pinned in the
[machine-readable receipt](ci_owned_group_timeout_20261002.json).

Apple's public [XNU signal implementation](https://raw.githubusercontent.com/apple-oss-distributions/xnu/main/bsd/kern/kern_sig.c)
filters zombie members during group signaling and can return EPERM when no
signalable members remain. This supports a disappearing/zombie-group race as
an explanation, not a fingerprint of the CI kernel or proof of that runner's
exact process state. Persistent permission failures must not be ignored.

Add one small POSIX cleanup helper shared by the measurement and WGD wrappers:

- Signal only the group created by the caller's `start_new_session=True`.
- Preserve TERM, a bounded grace period, and KILL for remaining descendants.
- On a permission error, poll/reap the leader. Accept cleanup only when a
  fresh signal-zero probe raises `ProcessLookupError` for the absent group.
- Propagate permission errors if the leader is live, the group still exists,
  or the second probe is also denied. Leader exit alone is not group exit.
- Record the helper's source hash in new measurement/native execution reports.

The measurement function continues exposing `stop_owned_group` for existing
callers. The WGD grace remains one second; it may now finish cleanup earlier
when the group is proven absent. Cleanup wall time is not a scientific speedup.
No process census or host contention monitor is added to cleanup.

## Validation

The final panel comprises **55 passing local cases in 8.24s**, zero failures,
errors or skips: measurement 37, WGD eight, measurement replay ten. It covers
actual owned Python commands, preserved exit/log/no-overwrite behavior, the
existing native child that ignores TERM, synthetic signal failures at TERM,
probe and KILL, persistent denial, a live leader, and a surviving group after
leader exit. Reports assert the helper's actual source identity. Eighteen
new cases are included. Both standalone CLI help invocations exit zero.

For the same final 14 permission-error cases, bind only the two historical
functions from committed source 11d in an isolated Python process, without
editing files: six pass and eight fail. Failures include the four absent-group
cases and four expectations for a confirming probe before propagating denial.
This is a synthetic comparison, not a historical native rerun. Retain the
earlier 11-case pre-fix report (five pass/six fail) and 51-case post-fix report
as well; those use the earlier fixture before its additional branches. Counts
overlap and must not be added together.

Actual source-11d macOS Python 3.11 fast CI has 13,746 pass, 30 fail, zero
errors, 110 skip and 30 warnings in 386.22s. All previous 85 retained-record
cases pass there, including their wrapped publication-figure log line.
At 05:17:27 UTC, all five source-11d test jobs are terminal failures;
CPU-wheel/docs succeed. Only the Python 3.11 fast log is inspected. The new
cleanup has local validation only at writing; its macOS/full CI outcome is
not yet known. No run is manually restarted or resubmitted.

## Scope And Next Work

Historical results, failed receipts, archived wrappers, calibration bindings,
frozen scientific settings, plans, run identities/order and accuracy values
remain unchanged. No completed inference, bootstrap, scoring, rendering or
archive reproduction is repeated. Before any prospective execution uses these
changed wrappers, its reviewed source/runtime inventory must include the new
helper; do not repin historical receipts or silently reinterpret old timings.

This fixes a bounded cleanup race, not arbitrary escaped sessions, complete
process-tree accounting, macOS resource equivalence or all subprocess APIs.
Inconclusive zombie-only groups whose permission error persists still fail
closed. Other snapshot/raw-data/Linux-capability CI failures, executable
release portability, QfO uncertainty, source/data rights and public deposition
remain open. The current 55-file review archive is not regenerated or expanded.
Timing remains deferred, with no quiet-window question, host contention poll,
DGX access or action on unrelated jobs/services. The publication goal is active.

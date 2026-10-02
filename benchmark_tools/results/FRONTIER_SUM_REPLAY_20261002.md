# Exact Historical Frontier CPU Replay

## Reproduced Mismatch

The retained source-c7b macOS Python 3.13 fast log rejects a portable
archived CPU replay. Reproduce the cause using the installed system Python
3.12.3 without installing or upgrading anything. Across the three retained
method measurements, sixteen fields differ: only outside-frontier CPU and
root-minus-frontier residuals, with maximum absolute difference
1.7763568394002505e-15 seconds. All other screening fields and decisions match.
The same input reproduces exactly under local Python 3.10.13.

[Python 3.12 documentation](https://docs.python.org/3.12/whatsnew/3.12.html#other-language-changes)
records the change to compensated float summation. The helper uses builtin
`sum()` for two descriptive aggregates; its result can therefore change with
the interpreter even when the raw counters do not. This is observed rounding
drift, not new CPU usage, interference or method-performance evidence.

## Correction

Use explicit left-to-right addition in the existing scope iteration for
those two aggregates. Preserve the archived calculation, including integer
zero for an empty outside set. Do not add tolerance, replace reports, change
raw counters, alter scope inventory/order or weaken strict replay equality.
Per-scope CPU values, validation, screening decisions and admission flags
are unchanged. This intentionally preserves historical serialization rather
than claiming improved numerical accuracy. Other summations are untouched.

Add two synthetic aggregate regressions and three all-method exact retained
screen comparisons. The helper source changes; prospective execution recipes
must pin its new bytes after stabilization. Historical receipts, frozen
diagnostic scripts and private runtime copies must not be silently repinned.
This does not authorize reuse of a stale source-bound readiness inventory.

## Verification

**172 local tests pass in 13.70s**, zero failures/errors/skips, across frontier
validation, periodic/boundary replay, hierarchy, pressure and wrapper checks.
The earlier 59-case Python 3.10 pre-change pass overlaps; it is not a failed
control or another independent replication. The actual Python 3.12 before
receipt retains all sixteen mismatches. After the correction, complete nested
screening equality holds for all three methods under Python 3.12, without
tolerance or reinterpretation of any timing result.

A fresh standard-library Python 3.12 child also reproduces all three screens
from 15 copied files/2,079,752 bytes. Its original-checkout read canary is
rejected, subsequent original-path events are zero, all fourteen project
source origins are staged and child subprocesses are forbidden. Temporary
copies are removed. This Python-event guard is not OS containment or full
native/raw pipeline or cross-host restoration.

[Source/test/report/log pins and commands](frontier_sum_replay_20261002.json)
retain the before/after and copied evidence. Only one diagnostic helper and
two tests change relative to source-6645. Scientific implementation/settings/
scores, collectors, strict replay validator and archived input stay unchanged.

At 08:09:01 UTC source-6645 run 36981155399 has full/3.10/3.12 tests live,
3.11/3.13 queued, and Linux diagnostics/wheel/docs successful. This correction
is not in that source. Do not infer new macOS proof or complete CI, and do not
restart a handle. Timing remains deferred: no host-contention poll, scheduling
question, DGX access or unrelated process/service action occurs. Archived DGX
data are read locally only. Other-QfO uncertainty, rights, controlled resources,
full executable release and public deposition remain open.

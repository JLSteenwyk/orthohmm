# Offline Remote Command Paths

Offline scaling-plan construction must preserve target-host paths without
resolving them against the preparation host's filesystem. The shared native
OrthoHMM command builder resolves input/output/metrics/tree paths, which is
appropriate for its default local harness translation but can rewrite a remote
plan through local symlinks. A real temporary symlink reproduces the strict
scientific-command mismatch before this correction.

## Explicit Construction Mode

Add keyword-only `resolve_paths=True` to `native_orthohmm` and
`resolve_orthohmm_paths=True` to `configurations`. The defaults retain existing
local resolution, including symlinks. Only the offline remote planner passes
false. That mode preserves lexical target paths, requires absolute interpreter/
input/output/metrics/tree paths without parent traversal, and rejects non-boolean
modes. It does not inspect or authorize the remote machine. OrthoFinder
construction and scientific settings stay unchanged. No comparison guard is
relaxed.

The complete returned remote plan matches all **15 fields and 27 runs** of
the retained plan exactly. Its five additional on-disk source/input metadata
fields are not returned by `plan`; the regression checks that exact key-set
boundary, not arbitrary field stripping. Historical JSON and source receipts
are unchanged. The two preparation helpers have new source hashes; prospective
execution inventories must pin them after source stabilization. Existing frozen
helper copies are not silently repinned.

## Verification

The local pre-change symlink test fails with the same strict comparison error
seen in the prior macOS log. That log does not expose its differing path tokens,
so the symlink establishes a mechanism, not the original runner filesystem
trace. New-patch remote confirmation remains open.

Add 21 cases, including default local harness parity through symlinks, lexical
remote paths with an optional supplied tree, relative/traversing path rejection,
typed mode validation, the real symlink regression and full plan equality.
**120 local tests pass in 1.80s**, zero failures/errors/skips, across scaling,
Threadripper derivation and prior metadata modules. The intermediate 108-case
pass overlaps and is not another independent replication.

All **120 cases also pass in a fresh isolated Python 3.10 child in 2.05s**.
The copied tree has 1,829 files/11,960,087 bytes and 81 staged project module
origins. The original-checkout read canary is blocked, subsequent original
open events are zero, child subprocesses are forbidden and temporary staging
is removed. Shared installed dependencies remain external. This is Python
audit-event evidence, not OS containment, cross-host/native scientific
restoration, complete CI or controlled timing.

[Source, input, test, receipt and log pins](remote_native_paths_20261002.json)
retain the pre-change failure, overlapping intermediate pass and complete
local/copied checks. Scientific implementation, defaults, settings, scores,
retained input bytes and command payloads are unchanged relative to dabcbb2e.

## Preceding Remote Confirmation

At 08:49:18 UTC, source-dabcbb2e run 36985130604 has Python 3.13/full tests
failed, Python 3.10/3.11/3.12 tests live, and Linux diagnostics/wheel/docs
successful. Its actual Python 3.13 fast log has **13,912 passes, 10 failures,
118 skips, zero errors and 30 warnings in 405.24s**. All 69 cases in the five
previous metadata modules pass. This closes that bounded prior remote check,
not all 128 previous local cases remotely, full CI or the new path patch.
Observe the next automatic run after commit/push; do not restart handles.

Lineage raw replay still needs an explicit historical scheduler-path/source
snapshot. Its source guards remain intact; the new frontier helper is not a
substitute for archived deployed bytes. Other raw/provenance/platform failures,
other-QfO uncertainty, data/runtime rights, controlled resources, complete
executable release and public deposition remain open. Timing remains deferred
without a contention poll, scheduling question, DGX access, unrelated process/
service action or new biological computation. Publication readiness is not
established.

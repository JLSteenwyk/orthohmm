# Immutable Initial Observation And Explicit Abort Resolution

This follows the retained [pre-native failure](THREADRIPPER_SHARED_ATTEMPT_22413.md).
The original scheduler FAILED 1:0, failed preflight and missing native endpoints
are unchanged. OrthoFinder did not start; no retry, successful-native result or
resource value is inferred.

## Repair Scope

Change three execution helpers, not the frozen scientific method or inputs:

- `threadripper_environment_worker.py`: shared-host preflight retains an atomic,
  immutable `preflight_initial_process_sample.json`, rather than hashing the
  growing process stream. Repeated checks require that copy to match the first
  stream sample, revalidate the actual observer/parked native process and job
  membership, and reject already released/finished workers. Legitimate later
  samples, including an append in progress, do not invalidate the initial copy.
- `threadripper_panel_progress.py`: an explicitly resolved `not_started` outcome
  is excluded from comparative timing with unavailable native resources. It
  cannot masquerade as ordinary completion, authorize retry, enter the overhead
  panel or bypass the original failed scheduler/review decisions.
- `bind_threadripper_panel_history.py`: bind the actual pre-native audit, abort
  gate, failed release/preflight, unchanged original session and successful
  repair validation. Borrowed identities, invented resources and unvalidated
  repair claims remain rejected.

The periodic observer's pre-release start, 30-second target and 35-second bound
do not change. The whole process stream remains mandatory for post-run replay;
copying its initial sample does not validate later samples or certify isolation.
Capacity, boot/process identity, runtime and release freshness checks remain
required. Legacy exclusive preflight retains its original strict one-sample
behavior. No unrelated workload or service changes occur.

## Retained Validation

The first ten-module regression passes 576 cases in 26.75s. Add the full
shared-preflight test with a periodic append between its two observations,
then capture the final regression, JUnit counts and stable source/test pins.
The [validation receipt](threadripper_prenative_repair_validation_20261004.json)
is 7,660 bytes, SHA-256
`07d0b0bcaefc384f9ea5ba7e52083b09d263f341fe5888daf3d4645c11a57fc0`.
All **579 tests pass in 27.01s**, with zero errors, failures or skips.
Source files are checked before and after this captured regression. This is
software validation, not a production preflight or native timing measurement.

Tests cover complete/in-progress append, changed first sample/copy, exited
observer/native process, already released gate, failed repair, malformed/borrowed
identity, unsupported retry/admission and both original failure-resolution
contracts. The actual failed attempt is not rerun. The original stream pin still
fails after append; the new initial-sample pin remains stable.

## Real History Resolution

The [resolution source](resolve_shared_prenative_failure_20261004.py) first
checks the unchanged independently reviewed failure receipt and fresh `sacct`
corroboration of its retained terminal controller. It retains four category
decisions: runtime passed on recorded checks, environment failed, resources
unresolved because inference never ran, outputs-or-failure passed for the
documented abort. It then records an explicit no-retry resolution.

The [resolution](threadripper_prenative_resolution_22413.json) is 9,827 bytes,
SHA-256 `d35501e6a969043c494a4bfb20f8a993972d5aea14f899c1c764d553f5fdf399`.
The [resolved session](threadripper_prenative_resolved_session_22413.json) is
2,395 bytes, SHA-256
`dc463a3bf9f27c5035ffdd82dbacc8470fb6226dbc1daabbc05fc9c8dde20f45`.
Both are unchanged copies of their canonical work artifacts. Original metadata
and decisions are preserved, rather than rewriting an earlier session as passed.

The [real-prefix binding result](threadripper_prenative_resolution_result_20261004.json)
is 168,715 bytes, SHA-256
`b5a6217612e3d37dabb9fae5c0e5836efc47993d6d04b69d18a0a44e330f8c94`.
It binds all eighteen original identities and points to index 18, full
OrthoFinder 3.1.5 on four proteomes, repeat 2. Original index 0 retains its
post-native cadence exclusion; index 17 retains its pre-native exclusion with
no resource measurement. This establishes reviewed history position, not
native execution authorization. Supporting raw/work evidence remains local;
these receipts are not a portable standalone archive.

## Remaining Execution Work

The previous runtime inventory/source recipe still binds the old three helper
files. Refresh the affected current-source runtime lookup, source/resource
recipe and shared readiness prospectively before any new run. Extend the
continuation/reporting integration to consume the resolved abort without
requiring a fake successful-native `review_run_17/summary.json`. Current partial
tables and figure remain their dated seventeen-resource-attempt snapshot;
the abort must appear separately with null endpoints in subsequent reporting.
The old launcher is not sufficient for index 18 and must not be invoked blindly.

No timing job is launched here. The 27-attempt panel, final resource/manuscript
reconciliation, portable/versioned study release and full publication objective
remain incomplete. Shared-host contention remains annotated with unknown,
potentially method-dependent distortion; no quiet window or DGX is required.

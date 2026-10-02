# Boundary Collector Environmental Review

Fixed a concrete integration defect before launching the observer-control panel.
This is component/offline synthetic validation, not a real environmental
handoff, observer slowdown result or controlled-resource admission.

## Defect And Scope

The new engineering executor selected the right collector but passed both arms
to the same periodic pressure review. That reviewer rejects interpoint gaps
above the frozen pressure cadence limit, even under diagnostic-only PSI policy.
The boundary collector intentionally retains only two native points. A long
valid boundary command would therefore fail solely because of its intended
observation spacing. Short boundary fixtures did not establish this integration.

The executor now propagates the **validated engineering task arm** to post-run
review. Default/production requests still use the original periodic behavior.
Boundary review requires an explicit diagnostic-only native-pressure policy,
a matching boundary report/job/native clock/collector policy, and exactly two
contiguous raw points. It does not reinterpret an eligibility-based policy.

Both points still require valid reads, stable boot/group/CPU identities,
monotonic PSI counters, bounded point-read durations and native bracketing.
Only the inappropriate **interpoint periodic cadence requirement** is absent
for the explicit boundary arm. Whole-run host-process sampling, policy/image
preflight, foreign CPU limits, stream continuity/gaps, release safeguards and
configuration-byte checks remain unchanged. No collector or scientific work
is altered, and no additional PSI sampling is introduced into either arm.

Boundary pressure uses `threadripper_boundary_pressure_review_v1`;
combined review uses `threadripper_boundary_environment_review_v1`.
Both explicitly state `periodic_pressure_cadence_checked: false`.
Pressure output reports `interval_average_some_percent` and
`maximum_point_duration_s`, not a periodic maximum-rate/cadence claim.
The observed boundary gap remains visible. Native PSI magnitudes remain
diagnostics; parked preflight pressure limits are unchanged.

## Executed Checks

Source is frozen at `9643773a502413b82b566df4cd59cae5e437ec1f`.
Focused regression passes 285 cases; the 13-module wider panel passes
**643 in 40.13 seconds**, zero failures/errors/skips. Coverage includes:
long boundaries; all three PSI diagnostics; bad/missing/extra points; read,
boot, group, counter and duration errors; mixed roles/collector identity;
changed reports; outside CPU/processes; host sampling gaps; and both executor
arms' exact reviewer arguments. Existing periodic/v1/v2 behavior stays tested.

Execute the actual reviewer on three persisted synthetic JSON/JSONL cases
under the existing isolated test Python, not the native timing controller:

| Synthetic Case | Process Review | Combined Decision |
| --- | --- | --- |
| Periodic assumption, two points 100 seconds apart | Pass | Fail: pressure cadence |
| Explicit boundary, same native bytes | Pass | Pass within this fixture |
| Boundary with injected outside CPU use | Fail | Fail |

Each case has 51 process observations/50 intervals and two pressure points.
The nominal 98-second native span is **fabricated fixture data**, not an executed
inference duration. The injected outside CPU maximum is 0.5 core equivalent;
it is not a measurement of this host. The pass certifies neither actual
workload isolation nor native outputs, runtime or resource validity.

Separate stdlib readback checks all 37 unique pinned files, equal native
bytes across cases, sample counts, boundary/native spans, decision flags,
and independently recomputes CPU maxima 0/0/0.5. The first readback expected
exactly 0.3 seconds for a point; fixture float-to-integer construction yields
299,999,999 ns. Its assertion fails before receipt creation; correct only
the readback expectation and preserve all raw artifacts without replay.
The compact [validation receipt](threadripper_boundary_environment_validation_20261002.json)
is 5,733 bytes, SHA256
`13884f658c3bc2dcd05141dad15c490f7e497747dbcd093b07271617d12270d2`.

## Remaining Work

No scheduler job/allocation, native inference, live host observation,
environment worker or approved real policy/readiness is executed here.
Do not treat the synthetic pass as actual environmental handoff or complete
pair admission. Final runtime/source binding, real handoff, a quiet window,
all 54 prescribed engineering outcomes and separate 27 production runs
remain outstanding. Old startup inventories do not attest these changed
execution helpers. Original protocol/plan/calibration, scientific method,
scores, manuscript PDF and archives are unchanged.

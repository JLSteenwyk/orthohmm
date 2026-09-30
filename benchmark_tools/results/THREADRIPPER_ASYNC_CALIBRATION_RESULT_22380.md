# Asynchronous Host Observation Calibration

Job 22380, submitted once after prospective commit `cc028714`, completed
successfully: allocation/batch/native step all COMPLETED 0:0, allocation wall
44 seconds and native step 41 seconds. The submission requests no requeue.
The [terminal receipt](threadripper_async_calibration_terminal_22380.json)
and [independent raw-replay audit](threadripper_async_calibration_audit_22380.json)
retain actual results under the unchanged
[acceptance protocol](THREADRIPPER_ASYNC_CALIBRATION_PROTOCOL_20260930.md).

| Prespecified check or measurement | Actual result |
| --- | --- |
| Native exit | Zero, not timed out |
| Expected process/main thread identities | 160 |
| Complete interior samples | 30; none incomplete |
| Maximum interior gap | 1.000758503 seconds (limit 1.5) |
| Maximum point read cost | 0.094280583 seconds (limit 1.0) |
| Witnessed worker CPU | 853.260950 seconds |
| Native-task CPU bracket, including wrapper | 860.448829 seconds |
| Native CPU excess above worker witnesses | 7.187879 seconds; original bound passed |
| Simultaneously retained touched allocation | 2,147,483,648 bytes |
| Native-step lifetime peak, including launcher | 2,414,571,520 bytes |
| Native command wall | 31.661611597 seconds |

All seven independent checks pass. The audit replays retained raw points and
resource scopes rather than trusting the precomputed run summary. Afterward,
all 50 direct audit evidence pins were freshly rechecked and matched. The
original diagnostic failure 22378 and deferred overall failure 22379 remain
retained; neither is rewritten or upgraded by this new experiment.

The initial/final host inventories still bracket native execution, with three
snapshots. Periodic scans no longer block the one-second point loop. This
validates the identified cadence fix on the stated workload, not a bound on
causal workload slowdown. Host-process sampling is still observational,
non-atomic and liable to short-lived process races.
The three raw snapshots retain 3, 34 and 2 process-read errors respectively;
the summary's zero observation exceptions does not mean zero raw read races.

The host summary directly reports competing CPU, with a maximum observed
79.2804 foreign CPU-core equivalents. These are process-interval estimates,
not scheduler reservations or controlled timing. Exclusive Slurm did not
establish a quiet host. No unrelated process or service was stopped.

## Remaining Boundary

This is 32-process/four-thread enumeration, not 128-core scaling; Python threads
share a process GIL. Peak memory is an allocation lower-bound check, not exact
RSS or an overhead-subtracted measurement. All point/cadence/accounting limits
are unchanged. Passing them does not establish continuous containment, causal
observer overhead, the complete native environmental-review handoff, a reviewed
ordinary-service/process policy, source readiness or controlled comparative
timing. The 27 production identities remain unstarted. No OrthoHMM inference,
benchmark score or default setting changed. Full publication goal remains open.

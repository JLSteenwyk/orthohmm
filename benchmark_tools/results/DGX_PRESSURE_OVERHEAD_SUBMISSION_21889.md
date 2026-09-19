# Replacement Overhead Panel Submission

Pushed source/protocol17fb1b5 and verified deploymenta35ca99 before
submission. The full preflight exercised the actual frozen enumerator under
the exact pinned interpreter and checked all runtime, system and recipe
records plus original inputs. No inference outcome was inspected to change
the experiment.105focused tests and shell syntax validation passed.

Submitted fresh complete array21889 held at2026-09-19T04:35:40Z, with
eligibility04:37:40Z. It requests18tasks with concurrency1, exclusive
spark-7ff0,20CPUs/96GiB, one-hour limits and no requeue. All preparation
finished before submission, leaving at least two minutes before eligibility.
The original21869failure and occupied directories remain unchanged.

Recorder21890 runs onbizon with1CPU/128MiB,25-hour limit and24-hour capture
deadline, using the same committed capture source as the previous panel
(SHA-256 `c692cb0fc28084558fb4e4f2d21282a94a6a93ceb32255660bbb3e8d8d279c64`).
Confirmed RUNNING state and a successful first poll showing21889still held
before releasing at approximately2026-09-19T04:36:34Z.

- [First poll](dgx_pressure_overhead_first_poll_21889.json):
  `fe31890d1eeb8af191af7ea36a6643c89120a15fe3415aae3c56d35806738156`.
- Submitted shell source:
  `369c64e3e992014279c883478761b837c1b9c9e9da064057c7f718f266a024ce`.
- [Array allocation request](dgx_pressure_overhead_submission_21889.txt).
- [Live recorder allocation](dgx_pressure_overhead_recorder_21890.txt).

Scheduler text uses EDT, four hours behind the UTC times above. Recorder
files are under `benchmarks/work/pressure_overhead_scheduler_21889/`, with
job log `benchmarks/work/pressure_overhead_capture_21890.log`.

## Frozen Bindings

- Plan: `b644e165dbf4d0beabf1cf4d9b6c314de522e3ebd1b91598ebebea99094c8fff`.
- Recipe: `50fc3a14c4d5b53c4eb158a22efcfb8b289b346718fad3c78b300172e319d8a1`.
- Authorization: `17943df08dfe21bf1975ff0dcb29dc66a9507b156e4698b1672e7517cdfd85e6`.
- Runtime/core, inputs, native arguments, task order and overhead budgets
  remain unchanged from the [pressure protocol](DGX_PRESSURE_OVERHEAD_PROTOCOL_20260919.md).

DGX quiet window is active from release until local accounting confirms all
18tasks terminal. No SSH/SCP/remote log reads during that window. Retain all
failures and detailed terminal records, do not select reruns from partial
outcomes, and do not treat submission/native exit status as timing admission.

# Prospective Threadripper Pressure Roles

Status: 2 October 2026. Corrected prospective execution policy, not an approved
environment, live handoff, production attempt or controlled timing result.
Timing remains deferred. No quiet-window answer is needed for this work;
no current contention polling, DGX access or unrelated job/service changes.

## Selection Risk

The historical whole-run evaluator treated CPU, memory and I/O PSI threshold
exceedances as environmental failures. The combined process review propagated
those failures to the executor. That could exclude a method for its own
resource stalls, contrary to the [frozen distinction](THREADRIPPER_TIMING_AMENDMENT_20260928.md)
between native method behavior and outside interference. System PSI records aggregated task stalls;
its magnitude alone does not identify a competing workload. See the primary
[Linux PSI documentation](https://docs.kernel.org/accounting/psi.html).

This is a prospective correction before the 27 production runs, not a
retrospective promotion of shared-host timings. Old policies, protocol files,
manifests, calibration admissions and failed attempts remain unchanged.

## Explicit Policy Contract

New private execution requires both fields below in the otherwise fully
reviewed environmental policy. This fragment is not an executable policy:

```json
{
  "schema": "threadripper_environment_policy_v2",
  "native_pressure_role": "diagnostic_only"
}
```

All other identity, process/service, CPU and cadence fields remain required.
No real approved policy file or threshold values are supplied by this report.
Hybrids, unsupported versions and missing/wrong roles fail closed. Executor
and worker reject legacy policy on the private route before native work or
environment observation. The historical shared route and default retained
replay preserve v1 eligibility semantics; compatibility is not current approval.

During the native interval, valid pressure magnitudes and threshold exceedances
are diagnostics only. The v2 pressure review exposes separate evidence and
diagnostic verdicts; it does not reuse the old eligibility field with a new
meaning. Missing reads, changed boot/job identity, invalid counters, gaps and
unbracketed native execution still fail evidence review. Strict outside-process
policy and CPU checks, configuration/inventory identities and post-read hashes
remain admission gates. Before native launch, the parked-worker preflight still
enforces prospectively reviewed pressure limits; v2 responses disclose this.
No magnitude is clipped and no estimated overhead is subtracted.

## Retained Replay and Tests

Source commit `45d59897b5a8f8dc2ab9b0d83be1871f7e79c2e2` implements the
versioned roles. Commit `99ab40b7cc013f6647092990e332c4e852cc8d5f` adds the
private-route requirement. The final six-module panel passes 273 tests in
8.38s, with zero failures, errors or skips. Its 42,817-byte JUnit receipt is
`benchmarks/work/threadripper_pressure_roles_private_gate_tests_20261002.xml`,
SHA-256 `496d9243a660e3049d67a28d1485736eabecb9a9a4225f9799c355f9748e4ef2`.

Earlier panels pass 231 and 278 overlapping tests; do not sum them. An attempted
extension names a nonexistent test module and exits 4 without collecting tests;
retain that selector failure. The final count changes because legacy-private
combinations are replaced by valid route/policy combinations plus early
rejection cases, not because evidence checks are removed.

Read-only replay uses the retained 22380 calibration's 33 pressure points,
32 intervals and native timestamps. The legacy result matches the evaluator
from `0a3250f6adead11a3763bd9a56b2487506e346e6` exactly. With illustrative
zero PSI bounds (not scientific policy values), both roles retain CPU/memory/I/O
maxima of 10.660229%, 0.005199% and 0.004408%, respectively. Exceedance counts
are CPU 32, memory 2 and I/O 16. The v2 evidence verdict passes while its
diagnostic-threshold verdict fails; legacy eligibility fails. Maximum period
is 1.028547544s within the illustrative 1.5s cadence bound. In-memory missing
I/O reads, boot mismatch and removed-point mutations still fail. Raw inputs
are unchanged. This does not attribute the historical stalls to native work.

The separate readback validates 53 unchanged file identities and all seven
replay source records against their recorded Git snapshot. Two caller files
changed afterward to require v2 private policy; their current bytes and the
latest synthetic tests are separately bound to `99ab40b7`. Do not describe
the replay's older caller identities as current mutable workspace identities.
Pressure evaluation was not repeated for that readback.

| Evidence | Bytes | SHA-256 |
| --- | --- | --- |
| [Retained replay](threadripper_pressure_roles_replay_20261002.json) | 27,829 | `d4c5c2b71f174cd15f79b9062d17e27189a9cc2ef0a472f48659eb7751fb2809` |
| [Independent readback and private gate binding](threadripper_pressure_roles_validation_20261002.json) | 20,029 | `06b7a5d38f77db1e07304e5d3bf6ccc1d5710921b95deb42d508c307dcd1ade8` |

## Remaining Gates

No native workload, calibration, scientific benchmark or timing panel is
restarted. Passing pressure evidence alone is not a quiet-host certificate:
retained outside-work/read-race limitations are not cleared. Neither the
27 production identities nor the 54 observer-control tasks are admitted.
Full-scale causal observer overhead, final peak-memory/resource scope, a
reviewed real process/service policy, live private environmental handoff,
current source/runtime readiness and a verified quiet window remain required.
Old recipes and runtime/helper bindings cannot approve these changed helpers.
Scientific defaults, inputs, resource allocation and run order are unchanged.
No accuracy advantage, complete release or publication-readiness claim follows.

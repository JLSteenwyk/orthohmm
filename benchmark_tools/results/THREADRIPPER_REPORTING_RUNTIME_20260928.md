# Reporting Runtime Rebinding

The [v4 runtime binding](threadripper_runtime_binding_v4_20260928.json)
replaces v3 for the reporting-stage collector. All seven previous inventories
were freshly generated from their original explicit roots. Their entry sets
are unchanged. Exactly two records differ: `measure_threadripper_scaling.py`
and `replay_threadripper_scaling.py`, the intended collector/replay changes
from commit 3aa2393c. All other old records, including all 191 native baseline
paths, match exactly. A separate inventory adds `report_finalization.py`.
There are now 176,711 records. No scientific source, installed dependency,
frozen input or native command was changed.

Raw full inventories remain in
`benchmarks/work/threadripper_reporting_runtime_20260928`; they were not
constructed by substituting new file hashes into old inventories. The binding
pins those files, the unchanged native baseline, command plan and inventory
generator, and records the before/after changed entries. Explicit-tree
limitations remain: no hermetic OS snapshot, continuous mutation detection,
or inference-branch completeness is claimed.

Two fresh native Python lookup probes completed in
`benchmarks/work/threadripper_lookup_reporting_{a,b}_20260928`.
For both OrthoHMM and OrthoFinder, signatures match each other and the prior
accepted probe. File coverage remains complete for the observed imports and
mapped files; scientific package origins remain the frozen core and installed
OrthoFinder 3.1.5 root. The [v3 lookup receipt](threadripper_python_lookup_v3_20260928.json)
pins both probes and the new binding. The inspector itself is unchanged.

The fixture driver now pins that new lookup receipt. The driver is recorded
separately at submission and in its started record; changing its receipt
constant is not represented as a change to the native scientific code.
85 focused runtime, fixture, reporting, replay and wrapper tests passed before
the pin change; all 12 driver/runtime tests passed again afterward.

## Submitted Checks

Update: all three jobs subsequently completed and passed independent review;
see the [outcome record](THREADRIPPER_REPORTING_OUTCOMES_20260928.md). The
submission-time observations below are retained as history.

The [submission record](threadripper_reporting_fixtures_submission_20260928.json)
retains the exact commands, driver hash, binding and lookup references:

| Job | Method | Dependency |
| --- | --- | --- |
| 22367 | OrthoHMM high sensitivity | none |
| 22368 | OrthoHMM satellite_v2 | afterok:22367 |
| 22369 | OrthoFinder full | afterok:22368 |

Each is the existing 16-gene installation fixture with fresh persistent
output, tmpfs input, JIT and scratch directories. The fixed 26-hour allocation,
release-budget guard, pre/post full runtime and native lookup checks remain
in use. These are diagnostic jobs, not any of the 27 production identities.
At the first authoritative queue check, 22367 was RUNNING and the other two
were PENDING on dependencies. Do not infer completion from successful
submission or restart them because an observation times out.

Next: poll these exact handles to terminal state, independently replay v4
measurements (including reporting memory), validate native formats and compare
canonical outputs with 22363-22365. Preserve any failure and do not retry
automatically. Quiet-host eligibility, production orchestration, full-scale
observer validation and the publication requirements remain incomplete.

# Native Dual-Bracket Diagnostic Submission

Source/plan commit `0e59a3a` was deployed as a fresh recipe. Derived three
commands from pinned pressure-panel tasks 1, 3 and 8, changing only
run-specific output/cache paths and diagnostic task metadata. Native method
settings, input records, runtime manifests and enumerator remain unchanged.
The launcher recomputes the plan from its pinned baseline and protocol,
rejects any difference, and requires its collector/inputs in the recipe.
Forty-eight focused tests and batch shell syntax checks pass.

## Frozen Bindings

- Plan: `dgx_dual_native_plan_20260919.json`, SHA-256
  `029f19b0e21f356387ae4e2df50310a225cdba4f356037814646af57226e1113`.
- Recipe: [retained manifest](dgx_dual_native_recipe_20260919.json), SHA-256
  `bea6cbee0a0a22c70a26b14eff137c65d95cd6d885f857d570dd11e515a34b46`.
- Archive: `benchmarks/work/dual_native_recipe_v1.tar`, SHA-256
  `c042d4b27771087e1e21968fae109bb6e85be2046b7494af775c5e5fb991248c`.
- Batch script SHA-256:
  `2113f2ab99f96519edfa63ce71799d57d217eb7d22a674819fc49f0e9509220a`.
- [Protocol](DUAL_BRACKET_NATIVE_PROTOCOL_20260919.md) SHA-256:
  `488125e6e349a83965a6f8ef75dcc477ce4b8bba5f5a13c0dedb09e88c4e9460`.

Transferred archive hash matched; all 437 extracted files were checked
byte-for-byte. Remote selection preflight accepted all three exact tasks
and confirmed the new output root was absent. The launcher interpreter is
the same pinned DGX OrthoHMM Python environment used in the pressure panel.
Full runtime/input checks execute before and after each native command;
submission and selection preflight do not establish that they have passed.

## Jobs and Observation

| Job | Method | Dependency |
| --- | --- | --- |
| 21912 | OrthoHMM high sensitivity | Initially held, then released |
| 21913 | OrthoHMM satellite_v2 | afterany:21912 |
| 21914 | OrthoFinder full | afterany:21913 |

The DGX scheduler queue was empty before submission. All jobs request
20 CPUs/96 GiB, exclusive spark-7ff0, one hour, no requeue, native timeout
900 seconds. Submitted all three before releasing the first job; dependencies
enforce method order even after failures. No SSH/remote reads are allowed
from release until all three jobs are terminal. No partial native outcome
has been inspected to choose a repeat, setting or inclusion rule.

Controller-only recorder 21915 runs on bizon, 1 CPU/128 MiB, four-hour
capture deadline and 4:10 scheduler limit. Its six tests cover exact non-array
identity, malformed record rejection and successful complete capture.
Source SHA-256:
`acf957ad45198b94a1a528551caefa34d55b1a34ffd722a4f689a9fab280dbd1`.
It retains immutable polls and first detailed terminal records under
`benchmarks/work/dual_native_scheduler_21912/`; log is
`benchmarks/work/dual_native_capture_21915.log`.

Confirmed 21912 RUNNING on spark-7ff0 and 21915 RUNNING on bizon, with a
successful first controller poll. The other two native jobs remain
dependency-pending. No scientific timing or overhead outcome is admitted.
After all three terminate, collect/replay raw measurements and validate
native artifacts, output equivalence, recipe/input/runtime and scheduler
identity before interpreting any diagnostic outcome.

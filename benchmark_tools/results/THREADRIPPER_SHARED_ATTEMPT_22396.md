# First Shared-Host Attempt: Native Success, Monitoring Failure

The frozen index-0 run, OrthoHMM high sensitivity on four proteomes, repeat 0,
ran as job 22396 on the local Threadripper. Native inference exited zero;
the enclosing job exited 1 after its post-run process-stream review failed.
These are distinct outcomes, not a failed orthology inference or an isolated
performance measurement. Preserve both, without rerunning or admitting a
comparative timing from this attempt under the original protocol.

## Independently Checked Evidence

See [machine-readable outcome](threadripper_shared_attempt_22396.json),
[review source](review_shared_threadripper_attempt_20261003.py) and the original
[live launch receipt](threadripper_shared_launch_22396.json). Raw evidence stays
in the retained work/results directories referenced by those records.

| Endpoint | Retained Observation |
| --- | ---: |
| Input proteins | 73,266 |
| Output orthogroups | 35,242 |
| Native wall seconds | 369.748703027 |
| Native CPU seconds | 9,315.218632 |
| Native-step lifetime peak bytes | 3,382,607,872 |
| Process snapshots / intervals | 14 / 13 |
| Process cadence violations | 1 |
| Maximum process start-to-start interval | 43.166799229 seconds |
| Pressure points / intervals | 371 / 370 |
| Maximum pressure observation interval | 1.0222614 seconds |

Wall time is the native monotonic command interval. CPU is the uncorrected
native task-subtree counter bracket, including wrapper work. Memory is the
native-step lifetime kernel peak, including its launcher, not pure algorithm
RSS or a whole-job teardown peak. No background or observer overhead is
subtracted. These raw endpoints are verified; the original environmental
audit remains failed. No accuracy score is inferred from valid output counts.

Verify both pre/post runtime checks, the actual frozen argv/input universe,
native completion, resource replay and a complete native partition. Independently
replay the process and pressure streams: reproduce the original failure exactly,
not a passing verdict. Runtime, resource and output categories pass; environment
fails. The session is not eligible for automatic progression under the original
history contract.

Slurm's live controller later purges this terminal job. Retain its original
successful terminal query and corroborate it with a fresh accounting query:
parent and batch FAILED/1:0, native step 22396.0 COMPLETED/0:0. No job is restarted
because of controller retention or observation expiry.

## Cause And Repair

The collector takes the first host process snapshot, performs the environmental
release checks, and only then starts its periodic observer. In this run the
preflight delay adds about 13.2 seconds to the first 30-second sampling period.
The first gap is therefore 43.1668 seconds; all later periodic gaps are about
30 seconds and the final post-command gap is about 10 seconds. Scan durations
are 3.15-5.17 seconds. The failure is an observer scheduling defect, not rejection
of the approximately 54-56 observed competing CPU-core equivalents.

Fix the observer to accept a validated prior monotonic anchor, and start its
serialized background thread immediately after the initial host snapshot,
before release checks. Anchor its deadline to that snapshot's start. Keep the
original 30-second target and 35-second admissible bound; do not loosen either
after seeing this outcome. Tests reproduce the 13.2-second release delay and
verify the next scan is due at 30 seconds rather than 43.2 seconds. Future
anchors, restart, thread errors, stale release, denied release and unfinished
descendants remain rejected. Default unanchored observer callers retain their
interface; the paired boundary collector is covered by regression tests.

Seven-module regression initially has 357 passes and 42 failures: legacy
overhead fixtures assume the production tmpfs input has never been created.
Use synthetic, test-specific parent output/input paths while preserving native
work and retained metadata checks. Do not delete the real input or relax the
production fresh-path guard. The corrected panel has **399 passes in 26.31s**.

## Source Binding And Next Execution

The first audit source encounters a stale September-29 resource-protocol
collector pin. Preserve that failed review. Reconcile source pins to the
actually executed source recipe in a [separate binding](threadripper_resource_endpoints_shared_binding_20261003.json),
using [the retained binding source](bind_shared_resource_endpoints_20261003.py).
Primary/secondary resource scopes and failure handling do not change. Its
retrospective timing is explicit; it neither fabricates prelaunch provenance
nor changes the failed monitoring decision. The second review reuses the
already checked native audit and corroborated terminal record; no inference
or expensive runtime calibration is repeated.

The review completes before the two-helper scheduling repair. Original helper
bytes are retained at Git commit 0c833313 and in the executed recipe/runtime
inventories. Historical bindings are not current-source certificates after
the repair. Fresh current-source/runtime/resource bindings are required before
the next identity. Resolve the retained infrastructure failure explicitly in
panel history, with its original failed verdict and comparability limitation,
before advancing. Do not overwrite run 0, silently retry it, skip an unreviewed
attempt, change scientific parameters or disturb unrelated jobs.

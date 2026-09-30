# Prospective Native Boundary Fixture

One new integration fixture is justified by the newly implemented
[boundary collector and raw checker](THREADRIPPER_BOUNDARY_COMPONENT_20260930.md).
Do not repeat successful periodic calibration 22380 or scientific inference.
This protocol is committed before submission and permits one diagnostic attempt,
not engineering workload pairs or production timings. No automatic retry.

## Fixed Execution

- Host: local Threadripper `bizon`; no DGX access or service operations.
- Native command: the pinned local `/usr/bin/true` binary, with no arguments.
- Controller: existing private Python runtime; no shared-package installation.
- Helper inventory: all current top-level `benchmark_tools/*.py`, with exact
  bytes/checksums; pin the submission script and retained controller provenance.
- Native settings: 32 affinity CPUs (IDs 0-31), 64 task slots, 128 GiB, native
  timeout 85,800 seconds, one-second completion polling, common 30-second host
  observer, exactly two native points. No GPU resource requested.
- Scheduler: existing `gpu` partition on `bizon`, exclusive allocation, 64 CPUs
  per task, 128 GiB, five-minute allocation limit, no requeue, one submission.
- Fresh run directory: `benchmarks/work/threadripper_boundary_fixture_20260930`.
  Require original output absence; preserve all logs/receipts if anything fails.
- Release guard: explicitly absent. Do not fabricate a passing contamination
  review or claim that this tests the full environmental handoff.

The short allocation bounds this `/usr/bin/true` diagnostic. It does not provide
enough time for the frozen production timeout or validate scheduler-budget gates.
An exclusive allocation does not prove the host is free of non-Slurm work.
Contention is retained by the common host observer; it cannot invalidate an
accuracy claim because no scientific inference or accuracy endpoint is run here.
Do not stop, suspend or signal unrelated processes or change scheduler/services.

## Prespecified Checks

Before command launch, require the exact protocol digest, helper inventory,
binary/interpreter/script/provenance pins, local cwd/host and task resource
environment. Unset `LD_PRELOAD`, `LD_LIBRARY_PATH` and `LD_AUDIT` in the script.
After collection, independently replay actual raw boundary evidence with explicit
job, command, interpreter-invocation and worker bindings. Require native zero
exit, no timeout and exactly two native points. Reconstruct native CPU bracket
and step-lifetime peak using existing low-level readers, with no subtraction.
Recheck all raw evidence and prospective bindings after the audit.

The in-job component audit is provisional until separate terminal scheduler
readback confirms allocation, batch and native step completed 0:0 with no
restart/requeue and the frozen submission command/cwd/resources/time limit.
Retain the submitted job ID, scheduler stdout/stderr, command and source revision.
Run a fresh raw audit after the job is terminal, retaining both audits and
checking their equality and evidence pins. An observation failure is not a
terminal job; poll the original job rather than resubmitting.

On any native, collection, binding or terminal-audit failure, retain the attempt
and logs and stop. No reclassification, automatic retry, timing correction or
replacement of original evidence. All timing, workload isolation, scientific
output-validation and publication admission flags remain false on success too.

## Preparation Validation

226 focused tests pass, including 31 new driver tests. Synthetic tests exercise
driver-to-raw-checker composition, independent readback/no-overwrite, changed
source/policy/resource/interpreter bindings and retained failures with one call
only. An initial test-copy mistake and an interpreter-alias binding mismatch
were corrected before native submission. The driver now preserves the actual
interpreter invocation path while hashing its bytes; no checker guard relaxed.
Bash syntax and private-runtime module import/help checks pass.

Success proves only the specified new native component integration. It does not
validate full-method overhead, output equality between arms, continuous resource
containment, complete loaded dependencies, common monitor cost, environmental
release policy, a quiet host, the 54-task engineering panel, the 27 production
timing runs or publication readiness. Keep historical pending plan bindings
unchanged; any future execution recipe must prospectively bind this implementation
and its actual admitted integration result.

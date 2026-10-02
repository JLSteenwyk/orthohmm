# Cgroup And Tmpfs Test Boundaries

Synthetic worker orchestration must not read the host's Linux scope or counters.
Real tmpfs copy checks must report missing capability rather than silently use
a disk-backed replacement. This milestone changes only tests and documentation;
production collection, preparation and timing admission remain unchanged.

## Observed Evidence

The already retained source-3e126 macOS Python 3.11 full-test log shows nine
worker failures reading `/proc/self/cgroup` and eight preparation fixture errors
creating temporary directories under `/dev/shm`. The overall unit result remains
failed: 13,601 passed, 78 failed, 20 errors, 106 skipped and 30 warnings in
546.98 seconds. Reuse the [pinned log](ci_cgroup_tmpfs_fixtures_20261001.json);
do not infer sibling causes or repeat its download.

At 03:07:07 UTC on 2 October all five source-3e test jobs are terminal failure,
with wheel/docs success. Source-f41 has five live test jobs at 03:07:08 UTC.
Its Python 3.11 job API subsequently reports successful compiler, dependency,
checkout-kernel, HMMER and MCL setup steps at 03:12:21 UTC; tests remain live.
No source-f41 log or compiler/library versions have been inspected here, and
native regression outcomes must not be inferred from setup-step success.

## Scoped Changes

Nine worker cases explicitly request a fixture replacing only the cgroup-file
read and their host snapshots. Ordinary JSON/file reads stay real. The two short
subprocess cases still run and retain actual exit status, log and clock ordering,
but their resource evidence is explicitly synthetic, not a host observation.
Two new missing/inaccessible-cgroup cases verify that the unchanged worker
fails before readiness, observer gate, host snapshot or native launch.

The eight actual `/dev/shm` copy cases retain native filesystem operations and
skip explicitly where that directory is unavailable; no mount or fallback is
created. Resolve their temporary persistent root to preserve the existing
direct-path rule. Fifteen new portable layout cases exercise unchanged path
validation without creating the declared tmpfs path: all three method layouts,
outside/relative/wrong/traversing inputs, copy bindings, unknown method, wrong
native argv, tmpfs output and pre-phylogeny rejection. Layout acceptance is not
filesystem availability, memory charging or preparation authorization.

## Executed Validation

**59 focused Linux cases pass in 0.89 seconds**, zero failures/errors/skips,
including 17 new cases. Module counts are worker 28, Threadripper preparation 23
and related native preparation eight. All eight actual copy cases execute here.

A separate process deliberately removes cgroup-file and `/dev/shm` availability
before test-module import. The same scope gives **51 passes and exactly eight
native copy skips in 0.93 seconds**, no failures/errors. Worker/gate, unavailable
scope and layout guards still execute. This is Linux capability simulation, not
actual macOS execution, and the panels must not be added together.

Both JUnit reports and two test-source pins are retained in the receipt. The
three production modules remain byte-equal to the base commit. No workflow,
package, compiler, workspace library, scientific input/default/score or admission
rule changes occur. Tiny owned files/subprocesses are not a scientific benchmark
or controlled timing. No DGX, contention poll/question, unrelated job/service
action or shared-package upgrade occurs.

Commit/push and inspect actual new CI. New remote fixture confirmation, native
regression results, other utility/path/data failures, comparable resources,
remaining uncertainty, source/rights and the versioned release remain open.
The publication goal remains active; controlled timing stays deferred.

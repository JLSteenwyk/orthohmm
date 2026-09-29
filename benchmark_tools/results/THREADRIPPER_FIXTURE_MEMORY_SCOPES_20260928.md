# Diagnostic Memory Scopes

The [machine-readable summary](threadripper_fixture_memory_scopes_20260928.json)
extracts explicitly scoped peaks from the three retained collector-v5 fixture
replays. It rechecks the replay evidence identities, rejects incomplete gauges
and inconsistent scope/peak relationships, and leaves final whole-job peak
unknown. These small, shared-host fixtures are not comparative performance data.

| Diagnostic job | Native step peak (bytes) | Job peak before inference | Job peak after inference | Job peak through reporting | Final whole-job peak |
|---|---:|---:|---:|---:|---|
| 22373, OrthoHMM high sensitivity | 618377216 | 254279680 | 698740736 | 698740736 | Unavailable |
| 22374, OrthoHMM phylogeny | 617283584 | 254349312 | 697860096 | 697860096 | Unavailable |
| 22375, full OrthoFinder | 276701184 | 253775872 | 357433344 | 357433344 | Unavailable |

The native-step cgroup includes its launcher. The job cgroup also contains
preparation, controller and observer work. Each value is cumulative only up to
its recorded read boundary; the broader reporting read precedes subsequent
validation and teardown. Equal post-native/reporting peaks do not prove a final
peak. Neither peak subtraction nor addition identifies observer-only memory.
These are cgroup counters, not interchangeable with per-process RSS.

The successful native completion boundary helps delimit the native command,
but does not recover a missing final job read or missing Slurm usage record.
The earlier descriptor-polling experiment remains negative evidence and was
not rerun. No values were inferred from absent scheduler fields.

Nine focused tests pass, including invalid gauges, observation errors, wrong
parent/step scope and decreasing peaks. This is a reporting guard, not a new
resource-eligibility policy or a substitute for full-scale validation.

An asynchronous question requests approval to investigate and, if supported,
enable non-disruptive Slurm job accounting. No configuration, daemon, service
or unrelated workload has been changed. Without that approval, continue other
publication work while retaining the explicit final-accounting limitation;
do not manufacture a passing readiness record.

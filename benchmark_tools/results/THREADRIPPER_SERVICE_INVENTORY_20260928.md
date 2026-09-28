# Local Services And Timing Conflict

The new read-only service collector inventories loaded Threadripper service/
timer configuration and running service names without approving background
work. It reuses the existing fingerprint parser locally; it does not call the
DGX collector, open SSH connections, or require DGX permissions.

Two initial service observations, a fresh typed process/CPU observation, and
a final service observation are retained locally. The
[review receipt](threadripper_service_context_review_20260928.json) binds all
three raw files and reports the configuration comparison and process-group
CPU estimates. This is environmental evidence, not a scientific timing run.

## Configuration Evidence

| Loaded inventory | Units | Fragment/drop-in file observations |
| --- | ---: | ---: |
| System services and timers | 115 | 125 |
| User services and timers | 81 | 81 |

All three configuration inventories agree, with no observation errors and no
unit reporting a pending daemon reload. The 206 file observations resolve to
201 distinct files; aliases/templates can share files. Every final observed
file was rehashed after collection and matched its retained size and SHA-256.
All three raw evidence files were rehashed too. This removes a file-read
permission concern for these observed files, not every host configuration.

The collector requests only Id, LoadState, FragmentPath, DropInPaths and
NeedDaemonReload for loaded services/timers. It retains file identities, not
unit contents, ExecStart commands or environment values. External scripts,
EnvironmentFiles, shared libraries, manager runtime overrides, unloaded units
and transient scopes are not covered by these fingerprints. Service presence
or a readable unit file does not make a workload ordinary or permissible.
The observations are sequential, not an atomic configuration snapshot.

## Observed Outside Work

The fresh capture measured **92.571041 observed CPU-core equivalents**, with
zero sampling errors and three unmatched outside processes. Per-process CPU
intervals are not simultaneous and can miss short-lived work; this number is
not a scheduling reservation or a promise of future utilization. No native
job was running, so the diagnostic excludes only the capture observer PID.
The Slurm queue observations were empty.

41 user services were reported running, including 32 with fungal-analysis
names. The largest observed groups were:

| Service or scope | Observed processes | Average CPU-core equivalents |
| --- | ---: | ---: |
| fungal-baliphy-independent-chains-20260927.service | 17 | 15.9978 |
| fungal-selected-matched-intervals-20260928.service | 17 | 15.9929 |
| fungal-neocallimastix-guides-20260927-v2.service | 3 | 15.9827 |
| fungal-fastml-optimizer-refinement-20260928-v3.service | 9 | 8.0030 |
| fungal-full-polynomial-ml-20260927.service | 10 | 7.9965 |
| fungal-mafft-marker-trees-resumed-20260928.service | 5 | 7.9905 |
| IQ-TREE in a tmux scope | 3 | 7.6459 |
| fungal-domain-coordinates-20260928.service | 6 | 3.9851 |

The receipt contains exact cgroups, names and all groups, including the tmux
scope that has no service/timer fingerprint. Associations use the nearest
service/scope component of the retained cgroup path. Names suggest workload
purpose but do not classify every process scientifically. Substantial active
BAli-Phy, Python, IQ-TREE and gainLoss work is directly observed. None was
stopped, suspended, restarted, masked or designated approved background work.

The top three services account for about 47.97 observed cores. Reserving a
Slurm allocation or choosing different CPU IDs would not eliminate this
whole-host concurrency. The prospective quiet-window requirement still applies.

## Validation And Remaining Work

29 focused local collector and reused fingerprint-parser tests pass. Tests
cover the read-only command set, local-host restriction, command timeout and
failure retention, malformed configuration, unreadable files, configuration
additions/removals/changes, and missing inventories. Missing evidence cannot
produce policy approval. The old DGX module was not modified.

This inventory does not populate `ordinary_processes`, choose numerical
background bounds, create a readiness review, or change the executor's frozen
configuration references. A real policy must distinguish acceptable ordinary
services from active analysis, cover the relevant external configuration,
and be frozen before timing. Native handoff validation, whole-run policy
application, full-scale resource accounting and a verified quiet window remain
open. No new benchmark or native integration job was launched in this turn.

Collect another pair into a fresh output path:

```sh
python -B -m benchmark_tools.observe_threadripper_services \
  --output /fresh/threadripper-services.json
```

Do not replace the retained observation with a later cleaner sample. These
files document the conditions observed at their own timestamps; a new capture
is separate evidence, not a retry of an inference or a controlled timing repeat.

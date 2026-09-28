# Kernel Identity And Changing Display Names

The previous exact-name policy check exposed 47, 57 and 50 changed records
in fixtures 22367-22369. Fresh inspection of those retained snapshots shows
that all 154 changes affect only the name field, not PID, creation timestamp
or cgroup. These observations do not contain kernel-type evidence and are not
retroactively upgraded by this work.

Linux exposes `Kthread` in `/proc/PID/status`. The inspected
[Linux 6.8 proc implementation](https://raw.githubusercontent.com/torvalds/linux/v6.8/fs/proc/array.c)
derives it from PF_KTHREAD and obtains workqueue-worker display names through
`wq_worker_comm`. The [workqueue implementation](https://raw.githubusercontent.com/torvalds/linux/v6.8/kernel/workqueue.c)
includes the current workqueue/description in those labels. Thus a changed
display name alone is not proof of a new process or executable. These upstream
sources explain the mechanism; they are not an attestation of the host's exact
Ubuntu kernel build. The host reports `6.8.0-55-generic` and exposes the field.

## Capture And Policy Changes

The new read-only `observe_threadripper_process_identity` module extends the
existing observer's returned snapshots without modifying that observer or
the frozen collector. For each observed process it retains Pid, Tgid and
Kthread fields with observation bounds, checks creation time and cgroup
before/after, uses a fresh Process instance for final identity checking, and
rejects a changing user-process name. Missing or inconsistent fields and
disappearing processes are retained as errors, not inferred kernel identities.
No command lines, environments, full status payloads or signals are collected.
Existing capture paths cannot be overwritten; interrupted captures retain
their available observations and failure status.

Process policy v2 requires same-boot typed snapshots. An explicitly reviewed
`reviewed_kernel_thread` entry may change display name only when its PID,
creation time and cgroup match and both observations have Kthread=1. No
name-prefix inference or automatically approved kernel-process class exists.
Ordinary-background entries require Kthread=0 and still match names exactly.
Unknown processes, PID reuse, cgroup migration, counter regressions, missing
identity evidence and observation errors are not waived. Kernel CPU use stays
in the unchanged diagnostic; neither a policy match nor kernel identity
establishes absence of interference or timing eligibility. V1 policies retain
their strict-name requirement.

Process creation times and type are not executable or service-configuration
attestations. A user process retaining its name across exec remains outside
what this evidence establishes. Service/executable review is still required.

## Fresh Host Observation

The [receipt](threadripper_kernel_identity_observation_20260928.json) binds the
raw capture, three source files and three historical snapshot files, all
rehashed after use. The new capture took about 9.54 seconds wall time with a
three-second requested gap. Snapshot collection itself took 3.20 and 3.28
seconds. This is one busy-host observation, not a bound on preflight latency.

| Snapshot | CPU rows retained | Verified kernel | Verified user | Missing type | Total observation errors |
| --- | ---: | ---: | ---: | ---: | ---: |
| First | 2,034 | 1,745 | 288 | 1 | 3 |
| Second | 2,038 | 1,746 | 291 | 1 | 1 |

All 45 same-PID/creation-time/cgroup name changes in this capture have
Kthread=1 in both observations. Unlike the old fixtures, this is direct
same-observation type evidence. Missing-type records are disappearing USalign
and tectonic processes; two additional processes disappeared during the first
base snapshot. Errors were preserved, with no retry to obtain a clean sample.

The CPU diagnostic found 93.108860 observed foreign core equivalents and six
unmatched outside processes. It excludes only the capturing observer PID,
not an active native job; there was no benchmark allocation. An empty Slurm
queue was also observed in this turn, but does not establish isolation.
No quiet-host claim, approved background inventory or preflight pass follows.
No unrelated process or service was stopped or modified.

## Validation And Next Work

114 focused identity, policy, observer and executor tests pass. Controls cover
missing/duplicate status fields, kernel-like user names, actual kernel-name
changes, wrong-boot observations, reuse/migration, decreasing counters,
out-of-bounds evidence, missing type, failure retention and overwrite refusal.
CPU diagnostics remain visible even when the narrow name rule matches.

The live environmental worker, reviewed service/executable inventory,
deadline validation, full-scale observer/resource checks and final source
freeze remain unfinished. No execution permit or production timing was
created. The fresh host observation confirms ongoing competition, not a
reason to relax the quiet-window requirement. This change concerns the
future measurement workflow, not OrthoHMM's scientific settings or scores.

To collect a new read-only diagnostic into a fresh existing parent directory:

```sh
python -B -m benchmark_tools.observe_threadripper_process_identity \
  --output /fresh/typed-process-observation.json --interval 3
```

# Threadripper Accounting Support And Current Load

Read-only inspection on 29 September 2026 UTC (28 September local time).
[Machine-readable evidence](threadripper_accounting_support_20260929.json)
records installed plugin and manual hashes, selected live configuration,
version-specific upstream source and the retained local host observation.
No scheduler configuration or service was modified and no job was submitted.

## Current Feasibility

At 03:47:21 UTC, two process snapshots separated by a requested five seconds
measured 99.416759 persistent foreign CPU-core equivalents. The Slurm queue
was empty. Visible work included IQ-TREE, Python and BAli-Phy under non-Slurm
user services/scopes. No sampled identity or CPU-counter errors occurred.
This is a short read-only feasibility observation, not whole-run isolation
evidence or a production eligibility decision. All unrelated work was left
running. Requested coordination for a quiet window; no response yet.

The 1.48 MB raw observation remains local at
`benchmark_tools/results/threadripper_resumption_load_20260929.json`, with its
checksum in the support receipt. No command lines or environments were captured.

## Installed Accounting

Slurm reports version 24.05.2, `JobAcctGatherType=(null)`, collection frequency
30 seconds, `accounting_storage/slurmdbd`, `proctrack/linuxproc`, and
`task/cgroup,task/affinity`. cgroup configuration selects v2 and constrains
cores, devices and RAM. Both cgroup and Linux accounting plugins exist in the
configured `/usr/lib/slurm` directory. Presence is not proof of activation.

The installed manual documents task-level CPU/memory collection and permits
plugin changes for future steps while existing steps retain their old plugin.
Approval for a possible accounting configuration change remains pending.
No reconfigure/restart or privileged modification was attempted.

## Version-Specific Memory Limitation

Current [SchedMD documentation](https://slurm.schedmd.com/cgroup_v2.html#limitations)
describes peak-based memory accounting, but it must not be applied blindly to
this older installation. In the [24.05.2 source](https://github.com/SchedMD/slurm/blob/slurm-24-05-2-1/src/plugins/cgroup/v2/cgroup_v2.c#L2000),
the task accounting function reads `cpu.stat`, `memory.current` and
`memory.stat`; the file contains no `memory.peak` access. It assigns current
cgroup memory to the returned RSS field. The local `cgroup_v2.so` likewise
has `memory.current` strings and no `memory.peak` string. This is consistent
with that source behavior, not proof of exact source-to-binary identity.

**Inference:** enabling the installed collector alone is not a validated
solution for brief memory spikes or final whole-job peak accounting. Task
statistics must not be relabelled as a simultaneous job peak or summed across
steps. CPU accounting and terminal state could still be useful, subject to a
separate controlled validation after approval. Existing direct cgroup peak
reads retain only their documented observation scopes.

Next: retain the validated private runtime/fixtures, obtain workload
coordination, and resolve the terminal accounting observation boundary without
assuming that turning on this plugin supplies kernel peak capture. Do not
upgrade Slurm or alter unrelated services without explicit authorization.
The 27-run scientific panel and Threadripper-only decision remain unchanged.

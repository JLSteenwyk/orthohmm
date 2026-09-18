# Cgroup Frontier Engineering Probe

Purpose: establish whether disjoint scopes outside a target can be sampled
without double-counting nested cgroups. This adds observation capability;
it does not change the existing CPU timing thresholds or admit any run.

Read the target and every sibling along its ancestry to the cgroup root.
These scopes form a disjoint frontier. Enclose those reads with root CPU
counter reads, retain raw counters/timestamps, directory device/inode
identities and ancestor-direct process counts. Reject hierarchy changes,
overlapping scopes/windows, counter decreases, boot changes and malformed
paths. Preserve signed root-minus-frontier residuals. Ancestor-direct tasks
and transient scopes are not fully attributed by this design.

Initial check: one two-second read-only observation on spark-7ff0 targeting
`/system.slice/spark-7ff0_slurmstepd.scope`, with no scientific native run.
The SSH connection and observer remain active and can consume CPU outside
the target. Thus observed outside activity is expected, not a quiet-host
estimate or an explanation of prior intervals. Retain source hashes and
both snapshots; replay all arithmetic locally. If scope churn rejects the
sample, retain that failure rather than repeating until a clean sample.

No services or unrelated workloads are modified. No scientific benchmark
is launched or repeated. This first probe does not validate integration with
native timing, collector overhead, causal attribution, or non-CPU isolation.

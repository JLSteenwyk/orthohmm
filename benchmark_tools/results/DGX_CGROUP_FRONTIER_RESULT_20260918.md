# Cgroup Frontier Read-Only Result

Committed/pushed collector and protocol1ebd457 before one two-second DGX
observation. No scientific native job was launched. Target:
`/system.slice/spark-7ff0_slurmstepd.scope`. Both snapshots retained stable
scope device/inode identities. Raw counters and source hashes are in
[dgx_cgroup_frontier_observation_20260918.json](dgx_cgroup_frontier_observation_20260918.json).

The target accumulated0CPU-seconds. Disjoint outside scopes accumulated
0.014109CPU-seconds, including0.006625inuser.slice,0.002599intailscaled and
0.001759incontainerd. Root outer accounting increased0.02CPU-seconds,
leaving a signed0.005891root-minus-frontier residual. Scope reads are
non-atomic and differ in accounting granularity; these are observed counter
differences, not exact workload bounds.

The SSH session and collector were active. This verifies that the collector
can record real disjoint outside-scope activity and replay it locally; it
does not establish quiet-host performance, explain the earlier satellite_v2
flag, measure native observer overhead, or justify any timing correction.
Service CPU use here is not evidence those services caused past flags.
Ancestor-direct tasks and transient scopes can remain unattributed.

35focused tests pass, including live-evidence replay and source identity,
synthetic negative residuals, scope changes, inode replacement, decreasing
counters, read order, boot changes, symlinks and malformed paths. Further
integration requires contemporaneous native and host brackets plus bounded
controls; existing scientific timing admission remains unchanged.

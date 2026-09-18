# Counter-Based Native Pipeline Smoke

Freeze code and this protocol before execution. Execute exactly three
sequential exclusive Slurm tasks on spark-7ff0: OrthoHMM high-sensitivity,
OrthoHMM satellite_v2 with inferred phylogeny, and full OrthoFinder3.1.5.
Use the existing645-protein/eight-species missing20_20261101 fixture and
the original launcher-smoke specification, SHA256
8af8a049c69b7e479e8443e1e0e14b74c78d67b345c03d7f557768488d057fa2.
Only replace the output-root launcher_smoke_v1 with counter_native_smoke_v1.
Do not change scientific flags, input identities, input enumeration, native
tool paths or dependencies. Use20CPUs/task and96GiB, array concurrency1,
one-hour allocation limit, native command timeout900seconds, no requeue.

Retain the existing runtime/input verification and fresh-input preparation
outside the native command boundary. Replace process-tree sampling with a
one-second cumulative host/observer counter stream. Run each native command
inside a distinct20CPU srun step with its own before/after cgroup snapshots.
The native step shares allocation CPUs with the batch observer; no dedicated
observer core or overhead-free claim is made. Preserve GNU time output,
native logs, step logs, commands, raw counters and final controller status.
Snapshot reads are not atomic; do not subtract host/native counts to infer
foreign load. Memory is cgroup accounting, not RSS.

Keep all failures and do not retry selectively. Require successful native
exit and runtime/input checks before investigating output correctness;
exit alone does not validate the biological output. Separately inspect
full input membership and native orthogroup/pair/tree artifacts. Read errors
remain visible and cannot be ignored to claim controlled timing.

These are functional engineering smokes, not a scaling panel, comparative
speed measurement, calibrated observer-overhead experiment or independent
accuracy validation. All publication/controlled-timing admission flags remain
false. Do not reuse the historical monitor's overhead gates as evidence for
this counter collector. The previous27timings remain descriptive.

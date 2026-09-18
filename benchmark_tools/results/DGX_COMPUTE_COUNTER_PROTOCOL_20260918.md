# Fixed-Work Counter Probe

Engineering only. Before execution, commit this protocol and runner. Use
one exclusive spark-7ff0 allocation, two CPUs/task, 2 GiB, ten minutes,
no requeue. Execute six trials with monitoring flags off,on,on,off,off,on.
Each uses a fresh one-CPU srun step and four sequential child processes.
Every child SHA-256 hashes the same 16 MiB x-byte buffer 256 times, accumulating
the 256 digests into a final SHA-256 checksum. Work is fixed, not time-limited.

The batch observer polls every 0.2 seconds in both modes, taking host/own-scope
snapshots at each poll only when monitoring is on. Both modes retain one
pre-work observer snapshot and native pre/post snapshots. Retain observer
process CPU, native work wall time, child CPU/checksums, raw counter evidence,
source hashes, commands, read failures and terminal scheduler state.

Require all 24 children to complete with expected checksums in the native
cgroup, separate from the observer's batch step of the same job. Require no
counter read errors, native snapshots to bracket work, native cgroup CPU
increment at least 95% of summed child measured CPU, and native memory.peak
at least the known 16 MiB buffer size. These are coarse accounting positive
controls, not statistical accuracy bounds. Retain failures without selective
retry. Report medians and monitored/off wall change descriptively with no
acceptance threshold, inferred significance or general overhead claim.

This is not a native OrthoHMM/OrthoFinder workload, a multithreaded saturation
test, an exclusivity proof or a scientific timing run. Memory accounting is
cgroup memory rather than RSS. Do not subtract mismatched host/native windows
to manufacture foreign-load bounds. Keep publication_ready and
controlled_workload_verified false. A real pipeline smoke, prospective timing
inclusion plan and controlled scaling evidence remain separate requirements.

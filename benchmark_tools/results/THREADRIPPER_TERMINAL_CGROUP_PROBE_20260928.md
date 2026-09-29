# Terminal Cgroup Descriptor Diagnostic

Job 22371 was a single owned diagnostic, not one of the 27 timing identities.
It requested one CPU, 256 MiB and two minutes on bizon. The standard-library
worker published its actual Slurm membership, waited ten seconds, allocated
and touched 32 MiB, waited three seconds, released it and waited two seconds.
The observer ran outside the allocation and opened only read-only descriptors
for the exact job ancestor's events, peak/current memory and CPU counters.
No cgroup writes, counter resets, unrelated signals or scheduler changes.

The allocation and batch step both completed 0:0. A subsequent queue query
was empty. The observer collected 1,463 samples with a nominal 10 ms delay
between samples. There were **zero error-free empty-cgroup observations**.
At removal all four held descriptors returned errno 19 (`ENODEV`), and path
lookup returned errno 2 (`ENOENT`). Holding these descriptors did not preserve
readable terminal counters after removal on this host.

The first readable memory peak was 9,371,648 bytes; the last was 42,979,328.
The last readable CPU usage was 254,269 microseconds. These are observations
before removal, not verified complete-job final values. Do not use them as
new timing results or infer that nothing changed after the last readable
sample. The observer's own overhead is outside this diagnostic job scope.

[Result and source identities](threadripper_terminal_cgroup_probe_20260928.json)
bind the submission, worker identity, observations, logs and probe source.
The [compressed raw observations](threadripper_terminal_cgroup_observations_20260928.json.gz)
round-trip exactly. Six scope-validation tests pass. This was one execution,
without retry; the failed capture mechanism remains negative evidence.

The [kernel cgroup-v2 documentation](https://cdn.kernel.org/doc/html/latest/admin-guide/cgroup-v2.html)
describes `populated` and hierarchical peak counters. Their documented meaning
does not guarantee a userspace reader can observe the short empty-to-removal
window. Faster polling or repeated successful attempts would not establish
that guarantee either.

## Next Requirement

Do not integrate this polling probe as a production final-accounting gate.
A validated lifecycle mechanism must retain final counters before removal,
or an appropriately configured accounting collector must preserve them.
Either approach needs full-job scope, timing/CPU/memory checks and failure
handling; no scheduler modification is authorized or performed by this probe.
Existing native and reporting-window measurements retain their stated scope.
Quiet-host and remaining timing requirements remain separate and unfinished.

The probe source is retained for inspection; rerunning it submits a new
diagnostic job and is not needed to reproduce the saved observation analysis.

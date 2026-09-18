# Quiet-Control CPU Field Diagnostic

Replayed all three retained quiet controls from array21820 without changing
their counters, commands, thresholds or admission decisions. Input SHA-256:
`29d0268a44d77c23e4822a2a062cc9d1aa4fffe88279885e82aab0a32835441d`.

The satellite_v2 interval3 remains flagged. Its host outer window reports
3.53 user,11.36 system and0.02 softirq CPU-seconds; nice,irq and steal are
zero. The enclosing job-parent counter reports3.278004 user and11.264579
system CPU-seconds, totaling14.542583. Subtracting these fields yields:

| Diagnostic component | CPU-seconds |
|---|---:|
| Host user+nice minus job user | 0.251996 |
| Host system minus job system | 0.095421 |
| Host irq+softirq | 0.020000 |
| Total host minus job | 0.367417 |

The host read-window increment is0.01CPU-seconds. Neither the interrupt
fields alone nor that read-window increment accounts for the residual.
This narrows counter-level explanations but does not identify a process.
Host/cgroup accounting has different scopes, granularity and update delays;
even the user residual is not proof of unrelated user-space work. Negative
residuals in other intervals are retained. Overlapping windows must not be
summed. No timing correction or scientific admission follows.

A read-only journal query covering21:47:28 through21:47:45UTC returned only
eight Slurm launch/affinity messages (queried with a100-entry limit).
The selected user-manager/logind/Slurm query over the full array likewise
showed only Slurm entries. These observations are not a complete activity
inventory: a service can consume CPU without logging. No service was stopped
or reconfigured and no benchmark was rerun.

Result:[dgx_quiet_cpu_fields_20260918.json](dgx_quiet_cpu_fields_20260918.json),
SHA-256`76d748b5007a24f964d2273d2f5f7784bdc2f987c1ad537841d30325c54b031c`.
The result pins its source/helpers and replays the original screens before
decomposition.50focused tests pass, including exact retained-result replay,
changed-input/screen rejection, field identities and signed residuals.

Next attribution work requires contemporaneous counters outside the job
scope or a validated accounting explanation; retrospective service names
cannot supply it. The27scientific timings remain unadmitted.

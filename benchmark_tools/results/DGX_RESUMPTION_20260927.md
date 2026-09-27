# DGX Resumption Check

The user renewed permission to use the DGX on 27 September 2026 and explicitly
requested public internet searches for TreeFam/QfO files instead of contacting
anyone. This supersedes the earlier DGX deferral. No external contact was made.

Read-only SSH to `jlsteenwyk@10.10.10.2` succeeds. Host `spark-7ff0` reports
20 CPUs and approximately 119 GiB total memory, 115 GiB available, with
3.2 TiB free on its home filesystem. Slurm's spark partition is idle with
20 CPUs and 106188 MB configured memory. The existing
[replacement protocol](DGX_SCALING_REPLACEMENT_PROTOCOL_20260920.md) already
specifies exclusive 20-CPU/96-GiB allocations, not 32 CPUs/128 GiB. The earlier
conversation's larger requested resources must not be described as available
on this host. No scientific setting or allocation protocol was changed.

The first deployed observer was an older revision lacking unit-configuration
fingerprints. Its successful command exits were not taken as configuration
validation. Deployed the current `observe_dgx_environment.py` under a fresh
filename and ran its bounded read-only static snapshot. Its recorded source
SHA256 matches the local source. All observation commands returned zero;
the snapshot lists no scheduler jobs or containers and no loaded Ollama models.
These are point-in-time observations, not a whole-run isolation guarantee.

Retained current snapshot:
`benchmarks/work/dgx_environment_config_resumption_20260927.json`, 113034 bytes,
SHA256 `23539fb805ea7344982b806ea11caa05cf97988420976cdb341953322bb1eb3a`.
The remote copy is under the existing publication workspace as
`environment_resumption_config_20260927.json`. Raw host inventory stays outside
Git; no credentials or unit contents were printed or committed.

Exactly one unit-file fingerprint failed: PermissionError for
`/etc/systemd/system/dgx-dashboard-admin.service`. A direct readability check
also failed. `sudo -n -l` reports that a password is required; no password was
requested or supplied and no privilege bypass was attempted. Asked the user
for administrator-provided read-only access to that one file. No blanket
sudo access or service shutdown is needed for this visibility fix.

The Samwise user service is still activating/restarting. Its prior narrowly
authorized temporary stop/runtime-mask and restoration policy remains relevant
for a future timing session; no service was changed during this check. No
native timing run was submitted while required environmental evidence remains
missing. Deployment, policy and whole-run validation must precede controlled
timing admission. All historical results remain unchanged.

# DGX Held-Session Service Probe

## Scope

Executed the new `benchmark_tools/probe_dgx_service_session.sh` on
`jlsteenwyk@10.10.10.2` in a single SSH session per probe. This is a bounded
service-maintenance diagnostic, not a scientific timing run, whole-host
isolation certificate, or production launcher. The local Slurm `spark`
queue was empty before probing. No Slurm jobs were submitted or cancelled.

The script requires explicit opt-in identifying the user's September 20
approval, the exact host and UID, the known persistent unit hash, its
original enabled link, and absence of pre-existing runtime masks. It stops
only the approved service, adds the higher-priority runtime mask, and
checks mask/load/activity/PID state four times at five-second intervals.
Exit/TERM/INT/HUP cleanup removes only its own unchanged symlink, reloads
the user manager, verifies persistent configuration, and requests the prior
start behavior. It never installs software or changes the unit file.

## Executed Checks

`bash -n benchmark_tools/probe_dgx_service_session.sh` passed.

Normal invocation:

```bash
ssh -T -o BatchMode=yes -o ConnectTimeout=10 jlsteenwyk@10.10.10.2 \
  'ORTHOHMM_APPROVED_SAMWISE_STOP=20260920 bash -s' \
  < benchmark_tools/probe_dgx_service_session.sh
```

- Exit code 0.
- Before: 2026-09-20T14:06:34Z, loaded, activating/auto-restart, PID 0.
- At 14:06:34, 14:06:39, 14:06:44, and 14:06:49 UTC: masked,
  inactive/dead, PID 0; fragment was the `/run/user/1000/systemd/user.control/`
  mask. All script assertions passed.
- Cleanup: persistent fragment loaded, start requested, transient PID
  1810064. `probe_exit=0 cleanup_failed=0 prior_active_state=activating`.

Deliberate TERM interruption:

```bash
ssh -T -o BatchMode=yes -o ConnectTimeout=10 jlsteenwyk@10.10.10.2 \
  'ORTHOHMM_APPROVED_SAMWISE_STOP=20260920 timeout --signal=TERM --kill-after=45s 8s bash -s' \
  < benchmark_tools/probe_dgx_service_session.sh
```

- Outer timeout exit code 124, expected for this interruption test.
- Before: 2026-09-20T14:07:07Z, loaded, activating/auto-restart, PID 0.
- At 14:07:07 and 14:07:12 UTC: masked, inactive/dead, PID 0.
- TERM cleanup: persistent fragment loaded, start requested, transient PID
  1810343. `probe_exit=143 cleanup_failed=0 prior_active_state=activating`.

Both cleanup paths verified the original unit SHA-256:
`cd2e76207b56adc59dbdaba796995f59ea7dea1320abc97c8fd2ecf5015cb84a`,
the original enabled link, and removal of the owned mask.
These observations are transcribed from command output, not a retained raw
controller archive or a scientific measurement receipt.

## Remaining Work and Limitations

A separate follow-up connection confirmed the persistent unit loaded and
auto-restarting. At 14:07:20 UTC its journal reported
`ModuleNotFoundError: No module named 'scientific_openclaw'`. Thus directory
creation repaired CHDIR only; transient `active/running` immediately after
start did not establish application health. The trusted application source
is still needed, as requested from the user. No persistent mask remains.

The tested service sequence must still be integrated into the replacement
timing session together with controller-terminal verification, continuous
environmental observation, and the frozen policy/preflight receipts. In
particular, interruption of `sbatch --wait` is not termination of the
underlying Slurm job: do not transplant this probe's unconditional exit
cleanup around an asynchronous job and assume isolation persists. A
production session must track that same job to terminal state or explicitly
record environmental interruption without resubmission or admission.

This short probe does not prove suppression throughout a 24-hour run,
absence of other workloads, performance neutrality of observation,
recovery after SIGKILL/host failure, or behavior when another actor changes
the unit/mask concurrently. It changes no frozen recipes, scores, or timing
admission decisions. Full replacement scaling remains unexecuted.

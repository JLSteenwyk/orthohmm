# DGX Samwise Service Maintenance

On 2026-09-20 the user requested repair of the missing working directory and
explicitly approved temporary stopping/runtime masking of this one user
service for dedicated timing, followed by restoration. Host:
`jlsteenwyk@10.10.10.2` (`spark-7ff0`). No other service was changed.

## Findings and Changes

- The enabled unit failed with `200/CHDIR` because
  `/home/jlsteenwyk/Desktop/Samwise` was missing. Created that directory.
- `/usr/bin/python3.12` returned `None` for
  `importlib.util.find_spec("scientific_openclaw")`. Directory creation alone
  does not repair the application. Bounded home-directory searches found no
  application source; this does not establish absence everywhere on the host.
  The user has been asked for the trusted installation or repository.
- Stopped `samwise-daemon-samwise.service`. The usual
  `systemctl --user mask --runtime` created a mask under
  `/run/user/1000/systemd/user/`, but the persistent home-directory unit
  took precedence: inspection still showed `LoadState=loaded`.
- Added a `/dev/null` symlink at
  `/run/user/1000/systemd/user.control/samwise-daemon-samwise.service`,
  which the inspected user-unit search path places above the persistent
  unit. After daemon reload, `LoadState=masked` and the fragment path
  resolved to this higher-priority runtime mask.
- The persistent unit was not edited, nor was its enabled symlink removed.
  Before/after SHA-256 of
  `/home/jlsteenwyk/.config/systemd/user/samwise-daemon-samwise.service`:
  `cd2e76207b56adc59dbdaba796995f59ea7dea1320abc97c8fd2ecf5015cb84a`.

These masks are temporary and do not survive removal of the user runtime
directory or reboot. Recheck actual load/active state before and during
timing. This maintenance is not evidence of whole-host isolation and does
not itself authorize admitting timing results. No replacement scaling run
was launched by this maintenance.

### Follow-up Verification: Mask Did Not Persist Between Connections

A later SSH connection, including a 15-second observation, found
`LoadState=loaded`, `ActiveState=activating`, `SubState=auto-restart`,
`MainPID=0`, and `NRestarts=2`. A subsequent check confirmed both runtime
mask paths were absent and `loginctl show-user` reported `Linger=no`.
Thus the earlier effective mask did NOT establish sustained suppression.
Stopped the service again, but do not claim it will remain stopped across
new user-manager sessions. No persistent mask or lingering setting was
introduced. The persistent enabled unit remains intact.

For timing, create both runtime masks within the single held SSH session
that spans the measurements, reload and verify `LoadState=masked` and
inactive state, observe throughout, and restore before closing that
session. Loss of the session/runtime directory invalidates the assumption
of suppression and must be handled as an environmental interruption, not
permission to silently resume timing. This lifecycle integration remains
unfinished. There are no surviving runtime masks at the last check.

## Required Restoration After Timing

Run as `jlsteenwyk` on the DGX, first verifying both paths are still the
recorded `/dev/null` symlinks. Stop if their identities have changed.

```bash
readlink /run/user/1000/systemd/user.control/samwise-daemon-samwise.service
readlink /run/user/1000/systemd/user/samwise-daemon-samwise.service
unlink /run/user/1000/systemd/user.control/samwise-daemon-samwise.service
systemctl --user unmask --runtime samwise-daemon-samwise.service
systemctl --user daemon-reload
systemctl --user start samwise-daemon-samwise.service
systemctl --user status samwise-daemon-samwise.service --no-pager
```

If the masks have already disappeared with the runtime directory, do not
run the unlink command against a nonexistent path; inspect the current
unit instead. The original configuration remains enabled. Restoring its start behavior
without installing the missing application will likely restore a failure
loop; do not describe that as a successful application repair. Complete
the trusted-source repair before claiming the daemon works. Keep the new
working directory, as explicitly requested. Restoration has not yet been
performed, and no application package has been installed.

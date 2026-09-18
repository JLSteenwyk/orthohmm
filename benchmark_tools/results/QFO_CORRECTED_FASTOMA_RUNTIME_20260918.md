# Corrected FastOMA Runtime Identity

`inspect_fastoma_runtime.py` completed two read-only inspections on bizon.
The reports are byte-identical. No container was launched, no package was
installed, and no Docker configuration or unrelated job was changed.

The retained report is `qfo_corrected_fastoma_runtime_20260918.json`, SHA-256
`c194ac79775f17d1067bdb5fefdb406e74edab6a61c976f8203cdb04023dbb99`.
The fresh recheck is
`benchmarks/work/qfo_corrected_fastoma_runtime_recheck_20260918.json`.

## Observed Identity

- 696 inventory entries across the installed Java tree and Nextflow capsule.
  No symlinks lead outside those roots.
- Seven checksummed binaries: Docker client, dockerd, containerd, runc,
  docker-init, Bash and the Nextflow launcher.
- Exact Java build: `17.0.18-internal+0-adhoc.conda.src`, reporting the
  version date 2026-01-20. This is the installed Conda build, not a newly
  installed or upgraded Java runtime.
- Nextflow 22.10.8 build 5860, probed offline.
- Docker client/server 28.2.2; containerd 1.7.27; runc 1.2.5;
  docker-init 0.19.0.
- Docker default context, `runc` default runtime, systemd cgroup v2,
  overlay2 storage, Ubuntu 24.04.2 LTS, kernel 6.8.0-55-generic, x86-64.
- Immutable FastOMA image identity matches the earlier pinned asset manifest.

The first inspection attempt rejected the abbreviated Java version assumption
because the installed executable reports `17.0.18-internal`. It wrote no success
report. The exact installed build was then recorded and tested, with no runtime
changes. All 22 focused runtime/asset tests pass, including wrong-build,
remote-override and incompatible Docker configuration rejection.

## Scope And Remaining Work

These are selected installed-file fingerprints, not a hermetic host snapshot.
Host dynamic libraries and active daemon executable mappings are not fully
inventoried. Container image identity does not by itself prove historical
container identity. The separate actual resource probe supports per-container
limits; this inspection adds no aggregate memory or dedicated-timing claim.

Before corrected inference, bind this report to the fresh command, staged
corrected inputs and admitted tree, recheck identity, and keep launch caches,
logs and task outputs isolated. Repeat identity checks after inference. Native
completion, pair conversion and scoring still require independent admission.
No corrected FastOMA inference or accuracy result exists from this inspection.

# Fresh Offline Scientific Reader Environment

The [relocated reader check](RELOCATED_INDEPENDENT_READERS_20260927.md)
still used the main installation's site-packages. A fresh reader-only venv
now runs the same four independent validators using copied local wheels and
a [five-package hash lock](publication_reader_requirements_20260927.txt).
No original inference package or environment was changed.

The [receipt](fresh_reader_runtime_20260927.json) records the lock, copied
wheels, pip install report, installed package inventory, before/after payload
audits, source-export manifest, file trace and report comparisons. Installed:
Biopython 1.86, DendroPy 5.1.0, NumPy 2.2.6, pip 26.2.1 and setuptools 83.0.0.
Installation was offline, binary-only and hash-required, with no retries;
`pip check` found no broken requirements. All 2,700 audited package payload
files match their wheels before and after scientific readback.
Sixteen focused package-audit and source-export tests pass.

All four readers pass on the existing relocated 16-gene native output.
Scientific fields in their reports and final summary match the earlier
readback exactly; 136 referenced-file occurrences also pass checksum checks.
DendroPy changes from the earlier reader's 5.0.8 to 5.1.0, explicitly recorded
rather than hidden by the comparison. This is fixture-level agreement, not
proof of general version equivalence.

## Reproduction

Actual fresh root: `/tmp/orthohmm-reader-runtime-20260927`. The five copied
wheels live under `wheels/`; `requirements.txt` matches the committed lock.
Choose a new root and use a trusted pip controller with these commands:

```bash
/absolute/base/python3.10 -I -m venv --without-pip /fresh/reader/venv
/absolute/trusted/installer/python -I -m pip --isolated \
  --disable-pip-version-check --python /fresh/reader/venv/bin/python install \
  --no-index --require-hashes --only-binary=:all: --no-cache-dir \
  --find-links /fresh/reader/wheels --report /fresh/reader/install_report.json \
  -r /fresh/reader/requirements.txt
/fresh/reader/venv/bin/python -I -m pip check
```

The actual validation used the new venv's Python with `-I -B`, explicitly
adding only the exported 31-module reader root to `sys.path`. A clean
environment set `HOME=/tmp`, `PATH=/usr/bin:/bin`, `LANG=C.UTF-8` and both
OpenBLAS/OMP thread counts to one. `strace -f -s 4096 -e trace=%file`
followed file accesses. Neither the original project mount nor the main
installation's site-packages prefix appears in that trace.

An earlier dependency-discovery probe that prepended old site-packages before
the standard library failed on an obsolete `pathlib` backport importing
`collections.Sequence`. This terminal observation is retained; no backport
was deleted and no shared environment was modified. Fresh venv execution
does not use that workaround or the old site-packages directory.

## Limits

This is reader-only isolation, not a changed inference environment, rerun of
inference/scoring, full-dataset validation or cross-host reproduction. The base
Python and OS libraries remain shared. The package auditor retains its
generated metadata/bytecode and non-site-payload exclusions. Literal trace
prefix checks are not a sandbox or universal file-access proof. The fixture
contains no satellite merges. Runtime archival, redistribution review,
dedicated timing and other publication requirements remain open.

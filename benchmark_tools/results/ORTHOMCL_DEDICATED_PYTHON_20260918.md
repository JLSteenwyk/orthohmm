# Dedicated OrthoMCL Preparation Python

Created `benchmarks/work/orthomcl_python_env_20260918/` as a fresh venv with
no seed packages and no system site-packages. The base interpreter remains
the existing CPython 3.10.13; no existing interpreter or package environment
was modified. Installed only Biopython 1.86 and NumPy 2.2.6, matching the
versions used by the tested preparation helpers, from local downloaded
wheels with `--require-hashes --no-index --no-compile`.

`orthomcl_python_requirements.lock` pins the exact CPython 3.10 Linux x86_64
wheel hashes. These are platform-specific artifacts, not a universal lock.
The installer report records both distribution hashes and environment
metadata. `pip check` reports no broken requirements.

## Clean Startup and Runtime Inventory

`freeze_orthomcl_python_runtime.py` requires isolated Python, disabled
bytecode writes, an absent cache prefix, a venv with user-site disabled,
and exactly the two pinned packages. It imports the actual corrected BPO
wrapper and rejects third-party modules outside the dedicated site-packages
or outside the Bio/NumPy module namespaces. No unrelated editable-package
startup hooks were observed. The repository metadata can appear in a global
distribution scan after adding its source root, so the installed-package
check is explicitly scoped to the venv's site-packages directory.

Inventoried 2,928 entries covering the venv, base interpreter, base standard
library (excluding base site-packages), and mapped native files. No external
symlinks remain outside the inventory roots. Imported helper sources are
recorded separately. The inventory and all imported-source/mapped-file
records were verified unchanged after the native fixture.

This is not an OS-wide hermetic image. The shared base interpreter and
standard library remain dependencies to verify before and after production
execution. The frozen executor must bind its own helper sources.

## Native Fixture

Using the dedicated interpreter with `-I -B -X pycache_prefix=<absent-path>`,
ran `prepare_orthomcl_bpo_checkpoint.py` against the retained guarded parity
fixture. It produced six BPO records from eight HSP rows, with three query
ranges and seven offset entries including EOF. Native index construction
and complete validation passed. `cmp` confirmed byte identity with the
native BioPerl fixture BPO, not just matching record counts.

45 focused tests passed for package/root validation, runtime inventories
and corrected BPO admission binding. The real dedicated-environment native
fixture is retained separately from those unit tests.

## Retained Hashes

- Runtime report `orthomcl_python_runtime_20260918.json`:
  `3519dbf43a376c935cc76ab00559af045c3e94ce47bf421319b970e660c8cddc`.
- Fixture report `orthomcl_dedicated_python_checkpoint_20260918.json`:
  `e78ad409617a112d03a33e5faaa33329e6bbd91a59ee4f0d256062d1376ee5a8`.
- Installer report `orthomcl_python_install_20260918.json`:
  `71a3b442e7c4252894dfb5002775b262e92be915ca25f885002382c598dcfaf1`.
- Wheel lock `orthomcl_python_requirements.lock`:
  `6ed53531e8485ac3443124f264a3eb41c9857ae9513b67248cd768201be0adae`.

Next: enforce this runtime identity in the production wrapper, freeze the
executor and batch command, and submit behind corrected BLAST admission.
No corrected full-data BPO, inferred group or QfO score was produced here.

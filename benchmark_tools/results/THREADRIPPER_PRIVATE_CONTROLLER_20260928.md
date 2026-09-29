# Private Timing Controller

Installed a separate controller at
`benchmarks/work/threadripper_private_controller_20260928/venv` using the
reconstructed private Python base, not shared Conda. The
[installation receipt](threadripper_private_controller_20260928.json) and
[six-package hash lock](threadripper_controller_requirements_20260928.txt)
record the offline installation. The five patched reader wheels were reused,
plus psutil 7.2.2 from the retained timing candidate. Wheel hashes must match
the prior reader lock and psutil receipt before installation.

All four stages passed: venv creation, hash-required binary-only offline install,
pip check and imports of the timing executor, verified fixture driver, native
output validator and collector replay. The controller has Biopython 1.87,
DendroPy 5.1.0, NumPy 2.2.6, pip 26.2.1, setuptools 83.0.0 and psutil 7.2.2.
These reader dependencies do not replace native OrthoHMM's scientific packages.
The checkout's source-distribution metadata is reported separately, not counted
as an installed wheel. The declared import probe observes no module from the
shared Conda prefix and reports the reconstructed private base prefix.

An initial read-only probe using the existing patched reader failed on missing
psutil. No package was added to that reader or the inference environment.
The new controller resolves that concrete dependency; optional execution paths
may still require validation. Imports alone do not establish execution parity.

## Runtime Inventory

The [inventory receipt](threadripper_private_trees_20260928.json) binds a local
6,173,195-byte inventory covering 18,751 entries across the controller venv,
patched inference venv and reconstructed base. No external symlinks were found
in these three trees. The immediate [recheck](threadripper_private_trees_recheck_20260928.json)
matches exactly. The large raw inventory remains local at the recorded path;
it is not an archived release artifact.

The existing snapshotter excludes bytecode caches and Git directories. This
does not inventory OS loader dependencies, scientific source, OrthoFinder,
external executables, controller helper files or every import path. Those
scopes must still be included in the prospective binding. An absence of external
symlinks does not establish absence of dynamic-library dependencies.

Eleven focused tests pass across the controller selection guard and runtime
snapshotter. No existing environment changed, no inference or production timing
was launched, and no DGX access occurred. Next: combine this private inventory
with the remaining native/OS/helper scopes, establish repeated lookup evidence,
then run collector v5 fixtures. Final accounting and quiet-host checks remain
separate prerequisites. No new comprehensive security assessment is claimed.

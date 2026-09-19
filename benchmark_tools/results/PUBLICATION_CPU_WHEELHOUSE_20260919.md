# Hash-Locked CPU Installation

Extended the [CPU-wheel smoke evidence](PUBLICATION_CPU_WHEEL_20260919.md)
with a locally retained11-wheel installation set, totaling92,916,492bytes.
The wheelhouse is at `benchmarks/work/publication_cpu_wheel_f9b9ce0/wheelhouse`;
no binaries were committed, uploaded or cleared for redistribution.

Downloaded binary wheels for the ten exact dependency/bootstrap versions
from the earlier installation, and copied the unchanged local OrthoHMM CPU
wheel. Versions alone do not imply identical dependency artifacts: this new
wheelhouse has its own exact file hashes and platform tags.

- [Requirements with SHA-256 pins](publication_cpu_wheel_requirements_20260919.txt):
  `60adb55b3647bf785319aec1348486003a28876fb5b656983add37b4cc496e99`.
- [Wheelhouse inventory and install bindings](publication_cpu_wheelhouse_20260919.json):
  `21cb1e91a94da6c0a96f19bed8fd256bbacf30884199e0cfa5ff49be7f570c42`.
- [Second installed smoke result](publication_cpu_wheel_offline_verification_20260919.json):
  `b7fa7a030d812aa0f4d51131b32e65343ccc4b382db142eeaf6b620077db00b0`.

## Executed Installation

Created a second Python3.10.13venv with system-site-packages disabled.
Installed every requirement, including pip/setuptools, with forced reinstall,
package indexes disabled, no cache, binary-only inputs and mandatory hashes:

```sh
"$OFFLINE_VENV/bin/python" -I -m pip --isolated --disable-pip-version-check install \
  --no-index --find-links "$WHEELHOUSE" --require-hashes --only-binary=:all: \
  --force-reinstall --no-cache-dir --report "$ARTIFACT_ROOT/offline_install_report.json" \
  --log "$ARTIFACT_ROOT/offline_install.log" \
  -r benchmark_tools/results/publication_cpu_wheel_requirements_20260919.txt
"$OFFLINE_VENV/bin/python" -I -m pip --isolated --disable-pip-version-check check
```

Installation succeeded and `pip check` reported no broken requirements.
Every one of the11installer download URLs resolves to a distinct local
wheelhouse file; its reported SHA-256 matches freshly hashed bytes. The
inventory exactly covers the directory. The retained pip23report uses the
singular `archive_info.hash` field. Initial inventory-check attempts assumed
the newer plural schema and passed augmented metadata into an exact file-record
checker; both failed without producing an inventory. The successful check
uses the actual report schema and unmodified nested file records. No package
or original installer evidence changed to pass these checks.

The same verifier then checked installed package bytes, venv-only imports and
CPU library loading, and ran both modes from outside the source checkout.
Standard and high-sensitivity each again produce4groups covering38genes
exactly once, with unchanged input bytes and the same retained partition
hash as the earlier smoke. Both actual subprocesses completed successfully.
Raw install/report/venv/configuration and verification evidence remain under
`benchmarks/work/publication_cpu_wheel_f9b9ce0/` with hashes in the inventory.

## Scope

This is no-index installation evidence, not an OS-level network isolation
test, a hermetic compiler build or a universal lock for other platforms.
Python, OS libraries and compiler requirements are outside the wheelhouse.
The OrthoHMM artifact remains a host-optimized CPU development build, not the
frozen scientific executor. Its local version0.5.0must not be confused with
an arbitrary public wheel of that version; the hash identifies these bytes.
Rebuilding may produce a different wheel hash. Public redistribution review,
portable release packaging, full scientific reproduction and publication
readiness remain open. No running benchmark or DGX access was involved.

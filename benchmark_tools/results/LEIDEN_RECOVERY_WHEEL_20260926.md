# Verified Leiden Recovery Wheel

The retained Leiden 0.11 snapshot used for OrthoBench recovery can be sourced
from the CPython abi3 Linux x86_64 wheel distributed by PyPI. The downloaded
wheel is `leidenalg-0.11.0-cp38-abi3-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl`,
SHA256 `571a0934f831a69442d82889d319bdba93de924bd9e09b720cd8cbe6fdc08c17`.
Its size and digest match the acquired official PyPI release metadata.

The [acquisition receipt](leiden011_recovery_wheel_acquisition_20260926.json)
records the exact artifact URL, provider metadata, download command and wheel
identity. The [payload audit](leiden011_recovery_wheel_audit_20260926.json)
validates ZIP paths and the wheel's complete internal RECORD inventory,
SHA256 digests and sizes without importing or executing its code.

All 15 payload members match the validated private snapshot byte-for-byte:
Python/C++ sources, the Python extension, bundled libigraph/libleidenalg,
package metadata, wheel tags and license. The only excluded comparisons are
the explicit installer-generated `INSTALLER`, `REQUESTED` and `RECORD` files.
The latter is independently validated for the wheel, not asserted identical to
the installed RECORD. No extra or missing payload member is accepted.

Six focused tests cover valid archives, traversal, duplicate RECORD rows,
bad hashes, bad sizes and unrecorded payload. No shared or private package
environment was modified and no native inference was rerun.

The [single-component hash lock](leiden011_recovery_wheel_20260926.txt) can be
used to download this historical platform-specific artifact into a new directory:

```bash
python -m pip download --isolated --index-url https://pypi.org/simple \
  --only-binary=:all: --no-deps --require-hashes \
  -r benchmark_tools/results/leiden011_recovery_wheel_20260926.txt \
  --dest /absolute/fresh/leiden-wheel
```

This lock does not supply igraph, NumPy, OrthoHMM, native alignment/tree tools
or a complete publication environment. A clean installation and installed-file
audit, followed by runtime/inference validation, remain required. It is not
a recommendation that this historical dependency is current or secure, nor
evidence of portability to another platform or full redistribution clearance.

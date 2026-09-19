# Native Build Isolation And Portability Audit

Following the platform-tag fix, source inspection found that wheel builds
compiled directly into the source tree and copied every source `.so` into
the wheel. A missing compiler or skipped kernel could therefore leave an
older binary eligible for packaging. That is incompatible with trustworthy
build provenance even when the wheel has a correct platform tag.

`CustomBuildPy` now copies package sources first, removes inherited `.so`
files from its build output, and compiles there. Native binaries are no
longer a source-package-data glob. Existing source-tree binaries are left
untouched; editable/develop behavior is unchanged. This change does not
modify compiler flags, native algorithms, or running frozen executors.

Seven focused packaging/entrypoint tests pass. The new parameterized test
simulates both available and skipped compilation, checks that stale source
and prior-build binaries cannot survive in the build output, verifies both
builders receive the output directory, and preserves the source binary.
This simulation is not a complete compilerless runtime test.

An actual build, local venv reinstall and installed CLI inference completed.
The fixture's38genes/4groups match the preceding installed-wheel output
byte-for-byte. [Verification](publication_isolated_wheel_verification_20260918.json)
records wheel/setup/output hashes. Evidence remains under
`benchmarks/work/publication_package_e76a248/`: `wheels_isolated/`,
`pip-wheel-isolated.log`, `pip-install-isolated.log`,
`installed-inference-isolated.log`, `installed-verification-isolated.json`
and `native-dynamic-libraries.txt`. The venv inherits dependencies; this is
not a clean dependency-resolution test or a portable release.

## Remaining Native Release Requirements

`readelf -d` on the newly built libraries establishes that the three CPU
kernels depend on `libgomp.so.1` and `libc.so.6`; Viterbi also lists the
Linux loader. The CUDA library lists libc and the loader; this does not
exclude dynamic CUDA driver discovery or establish GPU compatibility.

The build still uses `-march=native` when accepted. `hmm_have_avx2()` reports
the compile-time `__AVX2__` macro, not a runtime CPU capability probe.
Loading a shared library successfully is not proof its instructions are
safe on another x86CPU. Baseline-ISA compilation or validated dispatch,
dynamic-library compatibility, and cross-platform tests remain release gates.
The pair-alignment loader also directly loads its native library; broad
compilerless fallback claims need end-to-end validation rather than inference
from the search engine's Numba fallback. No portable binary was published.

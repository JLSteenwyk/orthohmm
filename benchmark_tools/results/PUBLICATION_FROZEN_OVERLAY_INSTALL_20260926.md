# Frozen Scientific Source, Setup-Only Build Overlay

Prepared and installed a new local CPU artifact combining scientific sources
at `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806` with only `setup.py` from the
tested packaging revision `6fd6df19daba83ec6467b917988f99e27a95be14`.
This is not the original scientific executable or a public release. No frozen
checkout, original source archive, package version or scientific default changed.

The staging tool materializes the same 43 selected files as the frozen source
archive, preserves the historical setup separately and changes only setup.py.
Build and installation leave all staged source bytes unchanged. An independent
[installation audit](publication_frozen_overlay_install_20260926.json) compares
the staged tree and all 33 shipped scientific source members directly against
immutable Git blobs, and installed bytes against those same members.
Six experimental reference files are omitted by the existing package selection;
they remain in the source archive. All three native CPU libraries are included
and load successfully; no CUDA binary is included.

## Artifact and Installed Checks

- Wheel: `orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl`, 144,444 bytes.
- Wheel SHA-256: `cfdfde5ed1be29e4080dd3571f5c0fc5fc5f45b57ebe5ee9b3c7559096c3b93d`.
- Audit SHA-256: `467036ee88d5e6f15d322aea4c94397591ef7ba521ac9f8ddde4b7cb783b97d2`.
- Exact 11-package [local wheel lock](publication_frozen_overlay_requirements_20260926.txt).
- Standard and high-sensitivity each produce 4 groups covering all 38 fixture
  genes exactly once, with unchanged input bytes. Both partition SHA-256 values
  are `1115fd8193636510bbc8cc8462d1b874e2a50db662fc0a59d3552d811ffa0885`,
  matching the earlier development-wheel fixture, not proving general equivalence.
- `pip check` reports no broken requirements. All eleven pip-report URLs point
  to the dedicated local wheelhouse; every wheel hash is freshly checked and
  the wheelhouse inventory exactly matches installed artifacts.
- An isolated installed capability probe returns `hmm_have_avx2() == 0`.
  Its local `baseline_probe.json` SHA-256 is
  `7be799fe5e5077c383b01ad11317e95c7b7345ee3df3b88a8a8532b23187a722`.

Eighteen focused staging, scientific-content and wheel-audit tests pass. They
cover real Git staging, source changes, unexpected/missing wheel members,
remote or hash-mismatched installation artifacts, duplicate/extra wheels and
output preservation. The installed smoke checker is unchanged and retains its
generic development-wheel limitations: this new packaging artifact still is
not the historical benchmark runtime.

## Execution

Artifacts reside under `benchmarks/work/publication_frozen_overlay_20260926`.
Build used the existing pip 26.2.1 / setuptools 83.0.0 environment, no package
index and no build isolation. PATH excluded nvcc and the explicit baseline
target omitted host-native/AVX2 compiler flags. The local dependency wheelhouse
was retained from earlier installation work, not freshly fetched or claimed
universally secure/current.

```sh
/usr/bin/python3 -S -m benchmark_tools.prepare_frozen_build_overlay --repo . \
  --output benchmarks/work/publication_frozen_overlay_20260926
# From the generated source directory; BUILDER is an absolute Python path:
env PATH=/usr/bin:/bin ORTHOHMM_CPU_TARGET=baseline PIP_NO_INDEX=1 \
  "$BUILDER" -I -m pip --isolated --disable-pip-version-check wheel \
  --no-index --no-deps --no-build-isolation --wheel-dir ../wheels --log ../build.log .
```

The builder is `benchmarks/work/publication_baseline_patched_install/venv/bin/python`
under the repository root. Its pip was used with `--python` to install into a
fresh `venv_clean` created with `python -m venv --without-pip`. The successful
install used one dedicated wheelhouse, `--no-index --require-hashes
--only-binary=:all: --no-cache-dir`, the lock above, and separate
`install_report_v2.json` / `install_v2.log`. The first attempt's `install.log`
and original empty venv are retained: two find-links directories exposed an
older same-name/version OrthoHMM wheel and pip rejected its wrong hash. The
repair isolated the intended artifact, not weakened the hash requirement.

```sh
python -m benchmark_tools.verify_cpu_wheel_install --root . \
  --python benchmarks/work/publication_frozen_overlay_20260926/venv_clean/bin/python \
  --wheel benchmarks/work/publication_frozen_overlay_20260926/wheels/orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl \
  --output benchmarks/work/publication_frozen_overlay_20260926/smoke
/usr/bin/python3 -S -m benchmark_tools.audit_frozen_overlay_install --repo . \
  --directory benchmarks/work/publication_frozen_overlay_20260926 \
  --output benchmark_tools/results/publication_frozen_overlay_install_20260926.json
```

Use new output locations for a separately justified reproduction; original
evidence is not overwritten. Recompilation need not reproduce the exact wheel
hash; a different binary needs a separate lock and verification.

## Remaining Gates

This is same-host installed-source and small-fixture evidence, not frozen native
dependency equivalence, full phylogeny/benchmark execution, non-AVX2 hardware
testing or cross-platform portability. The build environment, compiler and OS
libraries are not a hermetic toolchain. No rights clearance, publication release,
external archive, controlled timing, high-CPM admission or biological accuracy
claim follows. No existing benchmark result was transferred to this wheel.

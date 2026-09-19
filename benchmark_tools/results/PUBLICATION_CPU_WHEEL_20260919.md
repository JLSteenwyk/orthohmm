# CPU-Only Wheel Installation Evidence

Built a new local artifact from development revision
`f9b9ce09b01d075f6d4f81d65d730511fb7d9bae`, using a clean detached worktree
at `benchmarks/work/publication_cpu_wheel_source_f9b9ce0`. This is not the
frozen scientific executor7f3a9e4, a versioned public release, or a silent
replacement of the earlier CUDA-containing wheel. No artifact was uploaded.

## Build and Contents

The existing setup was used unchanged. Build PATH was `/usr/bin:/bin`,
excluding `/usr/local/cuda/bin/nvcc`; its log confirms CUDA was skipped and
all three CPU libraries were built. CustomBuildPy removes inherited shared
libraries from the build tree before compilation. The source checkout's
tracked setup/package files remain unchanged after building.

Artifact: `benchmarks/work/publication_cpu_wheel_f9b9ce0/wheels/orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl`.
Size166880bytes; SHA-256
`9584ef2e2b4797e94f50f3b809377708c5aa502b542e6833ce3e58bc7a6c733f`.
All42ZIP entries have consistent RECORD coverage; the embedded project MIT
license matches the repository. The only shared libraries are
`hmm_viterbi.so`, `kmer_prefilter.so` and `pair_align.so`. Their audited
symbol tables lack the selected CUDA runtime witnesses; there is no CUDA
shared-library entry. Project CUDA source remains included as source.
This is not a complete static-code attribution or redistribution clearance.

Build command, run inside the clean source worktree:

```sh
env PATH=/usr/bin:/bin /home/bizon/anaconda3/bin/python -m pip wheel \
  --no-deps --no-build-isolation --wheel-dir "$ARTIFACT_ROOT/wheels" \
  --log "$ARTIFACT_ROOT/build.log" .
```

The host build backend is not hermetic. CPU flags include compiler-probed
`-march=native`/AVX2; the Linux x86_64 wheel tag does not establish compatibility
with all x86_64 CPUs, other architectures or OS versions. A portable release
policy and separate testing remain required.

## Fresh Environment Verification

Created a fresh Python3.10.13 venv with `include-system-site-packages=false`.
Installed the local wheel and explicitly pinned NumPy2.2.6, Numba0.67.0,
llvmlite0.49.0, python-igraph1.0.0, igraph1.0.0, leidenalg0.12.0,
texttable1.7.0 and DendroPy5.1.0. Pip used cached dependency wheels; their
versions and installation log are retained, not a complete dependency
wheelhouse or hash-locked resolution.

`benchmark_tools/verify_cpu_wheel_install.py` runs the installed interpreter
with `-I`, from a fresh directory outside the source checkout. It verifies
that the imported package and dependency distributions are inside the new
venv, checks every installed package file against the wheel bytes, loads
all three CPU libraries, then runs standard and high-sensitivity native CLI
fixtures. Both produce4groups covering all38input proteins exactly once.
Input bytes remain unchanged. Both partition hashes are
`1115fd8193636510bbc8cc8462d1b874e2a50db662fc0a59d3552d811ffa0885`.
Isolated Python ignores Python environment flags; no global determinism claim
is made from this fixture. External phylogeny tools/pipelines were not tested.

The verifier itself runs in the host test environment with Biopython; only
the package probe and CLI run in the clean installation. Twelve focused
partition/wheel-audit tests pass in0.19seconds. The actual two-mode installed
smoke test also completes successfully.

## Retained Evidence

[Machine-readable verification](publication_cpu_wheel_verification_20260919.json)
SHA-256 `9e6cd7205322a67168c0776c9deb1d4d8722fbc03f0aeaf54669a667c8375f2f`
includes wheel inventory/linkage, installed bytes, dependency versions, exact
CLI commands, inputs, outputs and logs. Files remain under
`benchmarks/work/publication_cpu_wheel_f9b9ce0/`:

- `build.log`: `798ded27827e6ce788647250f76798728b5db8e33217c019a33351ad9f4f11e6`.
- `install.log`: `89e4dee3e48ea29e279c8990fa1fab016c8fcf8a1b79eeb627640546220fcc3a`.
- `venv/pyvenv.cfg`: `304a0dd3c736510cfb83ee13bcaaaf6fc81ce202d890acdaf82d48c376308a6e`.

This establishes a bounded CPU-only installation path on this host. It does
not transfer scientific results to development source, prove biological
accuracy or portability, resolve all release rights, or complete publication
readiness. Frozen running benchmarks and the DGX quiet window are unchanged.

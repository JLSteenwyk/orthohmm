# Baseline CPU Build Verification

Development revision `6fd6df19daba83ec6467b917988f99e27a95be14` adds the
opt-in `ORTHOHMM_CPU_TARGET=baseline` build target. It omits `-march=native`
and `-mavx2`, compiling all three CPU libraries, including scalar Viterbi.
The default native build behavior and all frozen scientific executors remain
unchanged. Invalid target values fail explicitly before compiler probing.

Thirty-nine focused setup/publication-runtime tests pass in1.27seconds,
including six target/AVX2 flag combinations, four invalid target values,
and actual baseline compilation/loading. The initial new test used the
wrong capability-symbol spelling; it was corrected to the existing
`hmm_have_avx2`, not by changing the kernel. No scientific code changed.

## Installed Verification

Built from the clean detached worktree
`benchmarks/work/publication_baseline_wheel_source` with:

```sh
env PATH=/usr/bin:/bin ORTHOHMM_CPU_TARGET=baseline \
  /home/bizon/anaconda3/bin/python -m pip wheel --no-deps \
  --no-build-isolation --wheel-dir "$ARTIFACT_ROOT/wheels" \
  --log "$ARTIFACT_ROOT/build.log" .
```

Here ARTIFACT_ROOT is the absolute path to
`benchmarks/work/publication_baseline_wheel_6fd6df1`. PATH excludes nvcc;
the baseline setting itself does not disable optional CUDA compilation.
The build leaves tracked setup/package sources unchanged.

Wheel `orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl` is145187bytes,
SHA-256 `3381a678fea1e85709305555d4211821bf1f75fa1d7268708804295a682b0c19`.
A fresh Python3.10.13 venv installed it with no-index, binary-only,
no-cache options and eight exact pinned dependencies from the existing
local CPU wheelhouse. Pip check reports no broken requirements. This
installation did not use require-hashes; the prior native-wheel lock is
not valid for this different artifact. Install report/log remain local.

The unchanged `verify_cpu_wheel_install.py` verifies wheel contents,
installed package bytes, isolated import paths, native library loading,
unchanged fixture inputs and exact output coverage. Standard and
high-sensitivity modes each produce4groups/38genes. Both partition hashes
are `1115fd8193636510bbc8cc8462d1b874e2a50db662fc0a59d3552d811ffa0885`.
A separate isolated probe of the installed Viterbi library returns
`hmm_have_avx2() == 0`.

[Unchanged verification report](publication_baseline_wheel_verification_20260919.json)
SHA-256 `c6b95b5b0ee57680bfe2aa9b53061ca588ff8c63daa0814288f8aaa141ae1de1`.
Build log SHA-256:
`d85471ce40b995f591505974cc9088147ffc7618746b001f12fee6b464d5c9a5`.
Install log SHA-256:
`ee332077a70feb390cef645fae2f90334d60ca98532800d4ab6de04ff72151fa`.

## Limits

This demonstrates a same-host baseline-ISA build and installed smoke test,
not non-AVX2 hardware validation, cross-platform portability, a manylinux
build, numerical equivalence across all inputs, or full phylogeny inference.
Compiler defaults, OpenMP/system-library compatibility and dependency-wheel
requirements still constrain deployment. No binary was uploaded, no public
release was made, and no benchmark result or timing was transferred to this
development artifact. Full scientific/package completion remains open.

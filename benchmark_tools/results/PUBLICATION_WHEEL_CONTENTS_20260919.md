# Retained Wheel Contents and Release Review

Audited the retained clean-install wheel, not the obsolete artifacts in
`dist/`. No wheel was modified or uploaded.

- Path: `benchmarks/work/publication_clean_install_4fe2c28/wheels/orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl`
- Size: 499209 bytes; SHA-256:
  `b06264e52edb784a6524a1a3a97c2099dbbd2e3a52cb5cbe57ac32d56b059e5b`.
- All 43 ZIP entries are covered by RECORD. All 42 non-RECORD entries
  match their recorded SHA-256 and byte count; RECORD is correctly unhashed.
- The embedded project license matches the repository MIT license exactly.
- Four native libraries are included, along with their C/CUDA source files.
  No standalone competitor executable, Python dependency distribution,
  proteome or reference dataset is present in this inventory.
- NumPy, Numba, python-igraph, leidenalg and optional DendroPy are declared
  installation dependencies, not bundled distributions in this wheel.

## Native Runtime Finding

All three CPU libraries dynamically require `libgomp.so.1` and `libc.so.6`.
The HMM CPU library also lists the system loader directly. These libraries
are not separate entries in this wheel. This is not an inventory of an
eventual container, full environment or dependency wheelhouse.

The CUDA library defines local `cudaMalloc`, `cudaFree` and
`__cudaRegisterFatBinary` functions, while its ELF dynamic dependencies list
only libc and the loader. Together with the build command's lack of a
`--cudart` override, this is evidence of embedded static CUDA runtime code.
[NVIDIA's CUDA 12.9 compiler documentation](https://docs.nvidia.com/cuda/archive/12.9.0/cuda-compiler-driver-nvcc/index.html#cudart-none-shared-static-cudart)
specifies static runtime linking as the default. The currently installed
compiler reports 12.9.41; this observation alone does not prove the exact
compiler used for every historical artifact.

The only standalone license document in this wheel is OrthoHMM's MIT file.
Do not treat it as a license grant for embedded NVIDIA code. The
[CUDA 12.9 agreement](https://docs.nvidia.com/cuda/archive/12.9.0/eula/index.html)
lists the static runtime among distributable components, subject to the
agreement's conditions. Before public binary distribution, identify the
actual linked toolkit components/version and resolve applicable distribution
terms and notices. This audit does not conclude that distribution is
prohibited or that adding a notice alone establishes compliance.

This finding does not alter scientific results: the frozen publication
runtime is CPU-only. Keep the existing wheel as historical installation
evidence; do not silently replace it with a CPU build or call it cleared for
release. A future release artifact needs its own inventory and review.

A subsequently built [CPU-only development wheel](PUBLICATION_CPU_WHEEL_20260919.md)
has its own inventory and isolated installation evidence. It neither replaces
this historical artifact nor establishes a public release or general portability.

## Reproduction and Limits

```bash
python3 benchmark_tools/audit_wheel_contents.py \
  benchmarks/work/publication_clean_install_4fe2c28/wheels/orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl \
  --license LICENSE.md
python3 -m pytest -q tests/unit/test_audit_wheel_contents.py
```

[Machine-readable evidence](publication_wheel_contents_20260919.json) retains
every member hash, metadata requirements, full direct dynamic-linkage output
and selected runtime symbol witnesses. The audit uses readelf/nm without
loading or executing the native libraries. Four focused tests pass, covering
valid inventory, corrupt content, unrecorded content and license mismatch.
The repository `venv` lacked pytest; the successful test command used system
Python. No full-suite claim is made.

This is a byte/linkage inventory, not complete static dependency attribution,
a license compatibility determination, or clearance of benchmark containers
and datasets. Those broader release reviews remain open.

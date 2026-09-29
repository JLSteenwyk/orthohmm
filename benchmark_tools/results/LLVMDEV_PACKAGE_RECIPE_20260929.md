# Resolved LLVM Package Recipe

The [producing-wheel log](LLVMLITE_PRODUCING_BUILD_20260929.md) resolves
`numba/label/llvm_wheel::llvmdev` to `22.1.0 manylinux_1` on linux-64.
Downloaded that exact published package without installing or executing it.
The [inspection receipt](llvmdev_manylinux1_metadata_20260929.json) retains
the release-metadata identity, package identity, complete info-member inventory,
selected metadata text and comparison with the release-tag build script.

## Verified Findings

- Published size: 893,725,486 bytes; downloaded size matches.
- Published and observed SHA-256:
  `c8603c82c26fb6b65c7bace30eaef09ce7287770d8e187ca5c93ec58b53e8c2d`.
- Embedded index and rendered recipe agree on build number 1, `manylinux_1`.
- Rendered source is LLVM 22.1.0, SHA-256
  `25d2e2adc4356d758405dd885fcfd6447bce82a90eb78b6b87ce0934bd077173`.
- Rendered `patches` is null. The template contains only commented patch
  examples; it does not list the AArch64 Windows patch found at the wrapper
  release tag. Do not transplant that tag's patch list into this package recipe.
- Embedded Linux build script is byte-identical to the previously retained
  release-tag script. It includes two `sed` edits to `AddLLVM.cmake`, so a null
  patch list does not mean unmodified upstream source.
- The recipe retains resolved build/host dependency versions and configuration.
  Metadata inspection was bounded to 32 MiB decompressed, rejected duplicate
  or traversal member names, and used no archive-path filesystem extraction.

## Evidence Boundary

The producing build log identifies the package by name/version/build, not its
checksum. This acquisition verifies the currently published package against
its release metadata; it does not prove that the historic build consumed these
exact bytes or establish a reproducible LLVM or llvmlite build. The source
archive remains to be acquired and verified. Only the info archive was inspected;
the package payload and complete license inventory were not inspected. The
metadata's `NCSA` label is not a substitute for inspecting the source notices.

No scientific runtime, benchmark configuration, endpoint or result changed.
Raw package remains local under `benchmarks/work/llvmdev_manylinux1_20260929/`.
Next source-provenance action: acquire the exact source archive using the embedded
recipe hash, retain applicable notices and the actual source transformations.

# llvmlite Release Build Source Trace

The upstream `v0.49.0` annotated tag resolves to commit
`b5a0ba74ae0601806c0ac3964d746f6be5f6d7b4`. Retained API metadata includes
the tag reference, tag object and non-truncated recursive tree. Six selected
build/recipe/patch files were downloaded at that commit and checked against
their Git blob identities, with independent SHA256 records for retained bytes.
[Four-file receipt](llvmlite_release_build_20260929.json) and
[two-file supplement](llvmlite_release_build_supplement_20260929.json).

The [wheel LLVM recipe](https://github.com/numba/llvmlite/blob/b5a0ba74ae0601806c0ac3964d746f6be5f6d7b4/conda-recipes/llvmdev_for_wheel/meta.yaml)
selects LLVM 22.1.0 and source SHA256
`25d2e2adc4356d758405dd885fcfd6447bce82a90eb78b6b87ce0934bd077173`.
Its declared archive is the upstream release's
`llvm-project-22.1.0.src.tar.xz`. The recipe enables an AArch64 Windows
target-detection patch, also retained at the same wrapper commit. The build
script additionally edits two CMake conditions, so the archive plus named
patch alone is not a complete description of build-time source transformations.
No LLVM source archive was downloaded or built in this inspection.

The [wheel builder](https://github.com/numba/llvmlite/blob/b5a0ba74ae0601806c0ac3964d746f6be5f6d7b4/buildscripts/manylinux/build_llvmlite.sh)
can install supplied LLVM development artifacts or resolve major version 22
from the project's wheel-build channel. Consequently, the release recipe is
a source candidate, not an immutable dependency lock for the retained wheel.
The inspected workflow also accepts a build-run identifier and retains its
uploaded artifacts for seven days; no specific matching run/artifact digest
has yet been established. The observed static LLVM 22.1.0 binary configuration
is consistent with this candidate but does not authenticate correspondence.

## Validation And Limits

Only public HTTPS files were retrieved; no contact, installation or build.
All six downloaded source files match their recorded Git tree blob hashes.
An attempt to extend the already-created parent receipt correctly failed with
FileExistsError. The original receipt was preserved, the two already-downloaded
files were reverified without another download, and a separate supplement was
written. This failure is retained in that supplement.

Exact wheel-to-build linkage, resolved LLVM package bytes, compiler/toolchain,
complete patches and source transformations remain unproven. Do not label
this a reproducible-build attestation, complete corresponding-source bundle,
security finding or redistribution clearance. No benchmark or runtime changed.

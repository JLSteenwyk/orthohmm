# Trace From Publishing To The Wheel Build

The attested upload run's current artifact listing is empty, but its job
logs remain available. [Listing/job receipt](llvmlite_upload_trace_20260929.json)
and [retained log checksums](llvmlite_upload_logs_20260929.json) distinguish
these observations; missing artifacts are not treated as missing logs.

The Find Workflow Runs log identifies Linux x86_64 wheel run `31462634019`.
The Download Artifacts log records successful retrieval of that run's
`llvmlite-linux-64-py3.10` artifact. The separately listed Conda build run is
not substituted for this wheel-producing run. [Producing-run metadata](llvmlite_producing_run_20260929.json)
confirms release commit `b5a0ba74ae0601806c0ac3964d746f6be5f6d7b4`, with
successful CPython 3.10 build/validation/test jobs. Its artifact listing is
also currently empty.

The [CPython 3.10 build log](llvmlite_cp310_build_log_20260929.json), job
`93689064413`, records the actual channel fallback rather than a supplied
LLVM artifact. Conda resolved `llvmdev-22.1.0-manylinux_1` from
`numba/label/llvm_wheel/linux-64`. This is build **1**, whereas the separately
inspected release recipe declared build **0**. Consequently the release-tag
recipe cannot be assumed to supply the exact resolved package's recipe or
patch set. Next inspect that exact package's metadata and embedded recipe.

The log also records manylinux container digest
`sha256:0a42cb7e5f4ba6bbfb8d0a86d1aab0c8876ba9c3be16bd99360ae42bf010ec77`
and uploaded artifact digest
`d872723ce04f4305d04ff6a9157e69b88e01907f23f8c084ad1a8208a69491be`.
The latter identifies the artifact archive, not necessarily the inner wheel;
it must not be compared directly to the PyPI wheel hash as if they were the
same format. No expired artifact was downloaded or byte-compared.

## Scope

This supplies a public-log chain from the verified publishing identity to
a producing run and its named LLVM package. It does not cryptographically
bind the missing build artifact's inner wheel to the retained wheel, establish
all compiler/linker inputs, or complete corresponding-source clearance.
Read-only log retrieval used existing GitHub authentication without recording
credentials or signed redirect URLs. Redirect downloads received no bearer
token. Logs remain local with checksums in receipts. No package was installed,
build run triggered, runtime changed or benchmark repeated.

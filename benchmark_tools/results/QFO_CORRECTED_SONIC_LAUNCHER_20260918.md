# Corrected SonicParanoid Launcher

The guarded launcher binds the prepared command manifest SHA-256
`188efb1f18e8f0bd90967670d18015ad500a9c47480333b0b380dff7c10c4c59`,
resolver snapshot `b4ed57defaba417f6627b671df958a2d533f4438bb8ff35a0b95cb820fd61bc9`,
and both full runtime inventories identified in
QFO_CORRECTED_SONIC_PREPARED_20260918.md.

It verifies all corrected input/source records and exact directory
membership, full runtime trees, and effective bundled-tool resolution
before inference and again after native exit. Fresh input copies are
checked independently. Existing output roots cannot be resumed or
overwritten. Native success remains pending separate output admission;
no conversion or scoring is performed by this launcher.

The native environment has explicit PATH, C locale, hash seed 0, disabled
user-site and bytecode writes, an absent isolated bytecode prefix, and
single-thread BLAS/OpenMP settings. SonicParanoid retains its explicit
32-thread command. Loader injections, CONDA_PREFIX and arbitrary Python
startup paths are excluded. This preserves the biological-application
execution policy without claiming identical historical environment.

The first real preflight failed because the reused OrthoHMM environment
verifier expected unrelated user-site package visibility. The corrected
launcher now checks the same frozen corrected input manifests directly
and uses the SonicParanoid-specific full runtime inventories/resolver.
No input or SonicParanoid dependency requirement was relaxed, and no
native output was created by that failed preflight.

Validation: 25 focused tests passed, including changed/missing/extra input
files, modified staging manifest, dependency drift, native failure,
postflight failure, check-only behavior, allocation and overwrite gates.
Batch shell syntax passed. Real preflight outcome and submission are
recorded below only after observed completion.

The real isolated-environment preflight subsequently passed with status
`preflight_passed_no_inference`, exit zero. Native inference is authorized
from the committed pinned executor; score admission remains separate.

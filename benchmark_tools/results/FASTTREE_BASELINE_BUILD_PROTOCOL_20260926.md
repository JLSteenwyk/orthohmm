# Baseline FastTree Build Protocol

Prespecified before build execution. This is a separate runtime-portability
experiment, not a scientific default change or a replacement for the historical
FastTree binary or the binary currently used by job 22179.

Use the already acquired FastTree 2.2.0 source at upstream revision
`29c5e62fbcd93230ee325f9c6a17b81f00e3c72a`, SHA-256
`975202a6b74c9996af871404ff043bb2152edcbda539035662514bc12d1f3431`.
Preserve the source notice and separately pinned LICENSE. Version 2.2.0 uses
double precision by default unless USE_SINGLE is defined. Do not enable
USE_SINGLE, OpenMP, fast-math or native/AVX compilation flags.

Compile twice in distinct private directories, using the same relative source
name and `/usr/bin/gcc -Wall -O3 -finline-functions -funroll-loops
-march=x86-64 -mtune=generic -o FastTree FastTree.c -lm`. Each build has a
240-second limit. Record compiler/linker identities, exact arguments, minimal
environment, logs and both binary hashes. A difference is reported, not followed
by a flag search. Record help output, ELF metadata and runtime-library listing.
This does not by itself prove execution on older CPUs or a hermetic build.

Use build_a with the rebuilt MAFFT and installed frozen-source OrthoHMM wheel
on the unchanged 16-protein fixture. Compare the same five primary outputs to
the retained historical fixture without changing seeds or acceptance rules.
Byte differences are reported rather than treated as proof of scientific
failure or suppressed. Native failure remains failed without automatic retry.
Preserve all old installations, evidence and the running job. Do not commit
third-party source/binaries or claim redistribution clearance.

After successful execution, independently validate outputs with the existing
structure, sequence, event and hierarchy readers. Full-dataset equivalence,
performance, cross-platform execution and dependency/source obligations remain
separate requirements.

## Help-Validator Correction Before Fixture Execution

The two builds completed with identical hashes, but the driver failed before
inference because it expected `FastTree Version ...` in `-help`. The retained
source at line 1813 explicitly prints `FastTree ...` for that option, and the
recorded probe exited zero with the correct version and double precision.
Preserve that failed report. Correct the banner validator with regression tests;
do not change compiler flags, rebuild, or overwrite outputs. A separate verifier
will hash-check the retained binaries and execute the first fixture in a new
directory, then run all four independent readbacks. This is an explicit harness
correction, not an automatic retry of failed native inference.

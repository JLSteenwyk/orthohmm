# Retained llvmlite Build Fields

The [probe receipt](llvmlite_binary_build_probe_20260929.json) records the exact
isolated-interpreter command and observed fields from the already installed,
full-run-validated llvmlite binary. This is not execution of downloaded source.
The binary and `binding/config.py` and `binding/initfini.py` were checked
against their exact wheel members before loading, and rehashed afterward.
The wheel itself was hash-checked before and after. No runtime changed.

| Reported field | Value |
| --- | --- |
| llvmlite | 0.49.0 |
| LLVM version | 22.1.0 |
| LLVM linkage | Static |
| libstdc++ linkage | Dynamic |
| Package format | Wheel |
| LLVM assertions | On |
| SVML support | False |

The command uses `-I -B`, a minimal environment and a 60-second timeout;
it exited 0 with empty stderr. It queries existing configuration fields and
LLVM version, not inference or JIT throughput. No optional inspection package
was installed. The receipt preserves the command for independent repetition.

The acquired wrapper source's `ffi/config.cpp` implements the linkage,
format, assertions and SVML fields through compile-time definitions, consistent
with the inspected Python wrapper API. This bounds their meaning as reported
build properties. It does not independently reconstruct the build, authenticate
a compiler or prove complete transitive dependencies. Static LLVM components
need not appear as dynamic NEEDED entries in the prior ELF inventory.

This replaces an assumption from the source's default LLVM major with an
observation of the retained binary's reported version and configuration.
Do not equate the version tuple with an exact upstream revision, patch set,
source-to-binary attestation, security finding or redistribution clearance.
Next examine upstream build records for this exact wheel and acquire matching
LLVM source/patch candidates only with explicit provenance. Frozen scientific
scores, runtimes and timing requirements are unchanged.

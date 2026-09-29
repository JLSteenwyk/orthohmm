# Numba And llvmlite Source Candidates

The existing source-acquisition workflow downloaded the source archives
advertised for the exact retained inference wheels. Published wheel and
source SHA256/size pairs, root package metadata and before/after local wheel
identities all passed. No downloaded code was installed, imported, built or
executed. [Machine-readable acquisition](jit_source_candidates_20260929.json).

| Package | Source bytes | Regular files | Notice candidates |
| --- | ---: | ---: | ---: |
| Numba 0.67.0 | 2,836,515 | 969 | 2 |
| llvmlite 0.49.0 | 194,467 | 101 | 2 |

The archives and raw release metadata are local under
`benchmarks/work/jit_source_candidates_20260929`; only the inspection receipt
is committed. The acquisition used the official
[Numba release metadata](https://pypi.org/pypi/numba/0.67.0/json) and
[llvmlite release metadata](https://pypi.org/pypi/llvmlite/0.49.0/json).

Source SHA256 values:

- Numba: `cd75aa535b33fa05d9d930b1ae8af9f97a2881e96d72dfb38ec9b78284d9f851`
- llvmlite: `00f16db782f4a13c78c5804aedc434e46794a77e89999a168f9401106270e50a`

## Remaining Native Boundary

Both archives contain top-level and third-party notice files, whose hashes
and member paths are recorded. These are filename-selected notice candidates,
not a complete component-to-license mapping or redistribution clearance.

Inspection of the acquired llvmlite `ffi/CMakeLists.txt` shows an external
LLVM CMake dependency, default supported LLVM major 22, an override for the
version check, and static LLVM linkage by default with a dynamic option.
Those are source build rules, not evidence of the actual retained wheel's
compiler, exact LLVM revision, patch set or linkage choices. This wrapper
source acquisition does not supply a verified corresponding LLVM source
bundle or close the static-component provenance gap. Next review the retained
binary's build metadata and upstream build records before selecting that
external source; do not infer an exact revision from the default major alone.

## Validation

All 31 source-acquisition tests pass. Independently rechecked the input and
tool records plus both wheel, metadata and source records after acquisition.
Source-member inventories were computed without extraction or build-backend
execution. Existing release/scientific limitations and frozen runtime bytes
remain unchanged; this is acquisition evidence, not a security assessment.

```sh
python -B -m benchmark_tools.acquire_wheel_sources \
  --inventory benchmark_tools/results/integrated_wheel_elf_20260927.json \
  --package numba --package llvmlite \
  --output /fresh/jit-source-candidates \
  --receipt /fresh/jit-source-receipt.json
```

Use fresh output paths; preserve failed attempts rather than overwriting them.

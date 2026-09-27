# Selected Numba Component Mapping

Inspected the pinned Numba 0.67.0 wheel without loading or installing it.
The wheel SHA256 was checked before/after; the third-party notice length and
digest match the [retained inventory](integrated_inference_notices_20260927.json).
The [mapping receipt](numba_component_mapping_20260927.json) preserves exact
archive members, uncompressed lengths, hashes and notice line references.

| Declaration | Evidence in this wheel | Scope |
| --- | --- | --- |
| appdirs | `numba/misc/appdirs.py` | Vendored Python module |
| NetworkX | `numba/core/controlflow.py` | Algorithm attribution, not an installed NetworkX distribution |
| cloudpickle | Three `numba/cloudpickle/` Python modules | Vendored package files |
| CUDA half-precision headers | `numba/cuda/cuda_fp16.h` and `.hpp` | Bundled headers, not evidence of benchmark GPU use |
| pythoncapi_compat | `numba/pythoncapi_compat.h` | Bundled header, not proof of compiled incorporation |

The control-flow module itself names NetworkX and includes an upstream URL
with revision `858e7cb183541a78969fed0cbcd02346f5866c02`. This is retained
provider attribution; no upstream content comparison was performed.
The third-party notice refers to CUDA Toolkit 11.2.2 terms for its headers;
the header text also identifies NVIDIA. These observations do not determine
applicable permissions or prove the CPU pipeline loads a CUDA runtime.

The notice's `jquery.graphviz.svg` declaration describes `docs/dagmap/`
documentation assets. No wheel member contains that directory substring or
the `jquery.graphviz` name. This is an explicit filename-scope absence check,
not a proof that related code cannot occur elsewhere. A source-tree notice
cannot be treated as a one-to-one inventory of installed components.

Five selected declarations map to eight hashed members. CPython/Unicode,
NumPy and numba-cuda attribution remains outside this partial mapping, as do
exact upstream equivalence, native/static component coverage and build/source
correspondence. The review adds provenance evidence, not legal advice,
redistribution clearance, runtime changes or publication readiness.

Reinspection uses the receipt's pinned wheel and notice records: open the
archive with `zipfile.ZipFile`, verify the notice bytes, inspect the recorded
line headings, and hash the eight listed members without importing them.
The receipt records both exact NetworkX attribution lines and the searched
documentation substrings so their limited scope can be reproduced.

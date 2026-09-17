# Native Banding Fix

The development C source now applies the short-pair exception separately
to each SIMD lane, matching scalar C and JIT semantics. Previously a target
longer than50 residues in a batch could incorrectly enable narrow banding
for other pairs whose profile and target were both at most50 residues.
Scores could therefore depend on batch composition. The stale comment
claiming that the multipair driver ignored band_width was corrected too.

The change follows the isolated original-fixture and boundary/order
diagnostics in narrow_band_rescue_20260917.json and
narrow_band_boundary_20260917.json. It is a correctness fix, not threshold
tuning, a new accuracy claim, or a demonstrated runtime improvement.

`tests/unit/test_native_viterbi_banding.py` compiles current C source into
a temporary library and compares SIMD/scalar/JIT scores over78 pairs at
five widths, three input orders and1/4 threads (30 cases). The mixed-length
fixture remains mixed after the Python wrapper's target-length sorting,
and covers empty targets and partial SIMD batches. With unpatched source,
the first band1 case fails at16/78 scores; with the fix all30 cases pass.
The full unit suite passes1865 tests. Tests explicitly skip without Linux
GCC or an AVX2-capable native build; no unavailable backend is called passed.

The publication baseline is still pinned at7f3a9e40dd7e79f842cc2c11fb8b548f9a802806.
Its checkout has no tracked OrthoHMM changes, and its native HMM library
still matches the original SHA256
f21d902852d727d6270b65ba900db52e7216cd874ee2d2b89fa72b6f16dcd5c0.
No installed library, frozen method manifest, dataset or historical result
was replaced. The development fix requires rebuilding the native library
when packaging a release; that release is not yet prepared or validated.
Default-width64 scores were unchanged in both retained diagnostic fixtures;
this evidence is not generalized into arbitrary dataset equivalence.

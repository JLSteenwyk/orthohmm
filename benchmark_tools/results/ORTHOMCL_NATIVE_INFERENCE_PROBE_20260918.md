# OrthoMCL Native Inference Compatibility

## Scope and Outcome

The bundled small OrthoMCL 1.4 BPO/GG fixture was run through native mode 4
in three fresh directories using the guarded Perl launcher: untouched
serial source, the existing pair-parallel patch with one worker, and the
same patch with two workers. All three exit codes were zero, with identical
partitions (12 groups, 41 grouped proteins from 42 input proteins) and all
112 directed graph edges and their printed weights identical. Matrix index
files match byte-for-byte. Raw matrix bytes differ in target ordering;
the complete validated, canonical weighted graphs match.

This is fixture compatibility evidence, not biological accuracy, equivalence
on corrected QfO, or validation of every 64-worker schedule. Graph edges
are compared only to diagnose native computation, not scored as final
ortholog predictions.

All three partitions differ from the bundled precomputed output: its one
11-member group is split into groups of six and five in every current arm.
The complete memberships are retained in the report. Because the untouched
serial arm has the same discrepancy, this difference is not specific to the
parallel patch. Its historical cause remains unresolved; no parameter was
tuned to remove it. A case-insensitive scan of all three native logs found
no `error`, `warn`, `fail`, `uninitialized`, or `reusing` matches. That scan
is not proof of absence of all possible native diagnostics.

## Prepared Sources and Provenance

`prepare_qfo_corrected_orthomcl_native.py` creates isolated native source
copies, changes configured paths and thread count, and applies the existing
opt-in pair-parallel patch. Original sources are checked unchanged. The
corrected production copy is prepared with 180 threads and 64 planned pair
workers, but has not executed. Scientific defaults, including inflation
1.5, are unchanged. Fresh paths are required; implicit cache reuse is not
allowed by this preparation.

The existing runtime inventory is supplemented with `/usr/bin/mkdir`,
`/usr/bin/time`, `libselinux`, and `libpcre2-8`. The native mode-4 source uses
`mkdir` through a spaced `system (...)` call. These helpers are verified
alongside the existing Perl/MCL inventory before and after the fixture and
source preparation. This does not establish an OS-wide hermetic runtime.

Retained result SHA-256 values:

- `orthomcl_native_inference_probe_20260918.json`:
  `838eb833fe0b18b429db1c9bb1ad759313638dbbd2dcbcde7eb1c8270f21dc4e`.
- `qfo_corrected_orthomcl_native_sources_20260918.json`:
  `9ccd23429301d754f80fc72ae4ba4b28ef5036116e0670c5aebdb85a41b00ae6`.
- `qfo_corrected_orthomcl_system_helpers_20260918.json`:
  `79ef5d6297914cf43e318858de4de42e58f31d601159a9e27bb8515734efaa2c`.

Final native artifacts are under
`benchmarks/work/orthomcl_native_inference_probe_v3_20260918/`.
Earlier probe directories are preserved, including the first comparison
that required raw matrix byte equality and reported a difference.

## Verification and Remaining Work

The 25 focused tests cover isolated configuration, preservation of the
originals and scientific defaults, rejection of invalid resources and
shell-unsafe paths, patch application, native membership validation, and
weighted-graph comparison. Native execution is additional evidence rather
than a mocked unit-test claim.

Corrected BLAST must pass independent admission before full BPO conversion,
content audit, native indexes and inference. Production scheduler binding,
runtime checks, final-group admission, conversion and QfO scoring remain
required. Prepared sources and passing fixtures do not authorize treating
the corrected OrthoMCL result as complete or publication-ready.

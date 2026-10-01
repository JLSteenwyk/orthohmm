# Explicit Historical Exporter Binding

The failed partial-panel export is retained in
[its failure record](QFO_PARAMETER_PARTIAL_EXPORT_FAILURE_20261001.md).
This amendment changes visualization evidence binding only, not inference,
reference scores, bootstrap arithmetic, endpoints, multiplicity or defaults.

## Scope

The default export route still rejects changed source identities. A caller
must explicitly request `--historical-exporter-binding` for the known old
root `benchmark_tools/export_qfo_corrected_comparison.py` record, 8,652 bytes,
SHA256 `4239f8d263295b313c5cc27b866bf7c6494feac44a8fbfe18e160bd92ff69bfb`.
No other source or data record can use this route. Requesting the amendment
without that exact original record is an error.

Bind the original exporter to its direct retained historical copy under
`benchmarks/work/publication_qfo_cpm_candidates_admission_v3` and verify its
bytes against Git revision `5efb206b23a44f85386d8cb7e90f3e815a3d162a`.
Separately require the direct current exporter, 9,572 bytes, SHA256
`8fd211aeaaac5765b6f5ffa01812d8e61ac4134b8c1ceb1bada4f497df7d7aa0`,
and verify its bytes against Git revision
`01ac4f66b2ac4b3565e5617b7abfe789e7ffe9aa`.
The latter change adds recovered OrthoMCL reporting. Neither historical nor
current exporter is substituted for parameter scoring/bootstrap code.

Preserve the original declared input records in the new manifest. Record
effective checked records, historical copy, old/current Git revisions,
current source identity and source-binding implementation separately.
Deduplicate only equal effective identities; conflicting records fail.
Check every other input unchanged, including the exact uncertainty result,
numerical-reproduction receipt, native admissions, conversion, raw/reference
evidence and all frozen scientific helpers. Recheck bindings before native
table validation, before writing and after output generation. Do not edit
historical receipts or revert current sources. A fresh output directory is
required; failure must not trigger an automatic retry.

## Validation And Interpretation

Tests cover strict default behavior, missing/unknown originals, historical
and current file changes, both Git-blob mismatches, symlinks, scientific-data
and scientific-helper changes, conflicting identities and explicit metadata.
The real export still requires full input checking; passing fixture tests is
not scientific admission. Reuse the pinned independent numerical reproduction,
not a new bootstrap. All seven arm positions and 18 planned intervals remain,
with unavailable high-CPM explicitly null in the retained partial panel.

GO/EC/FAS are not F1; the mean is secondary. QfO is development-exposed,
FAS unseeded, and no new uncertainty, independence, equivalence, superiority,
default selection, controlled timing or publication-readiness claim follows.

# Simulation Native Runtime Correction

This execution amendment applies to both frozen simulation panels and is
defined before corrected-runtime inference or variable-length accuracy
inspection. It does not replace either panel's scientific analysis plan.

## Fixed Scientific Design

Retain source revision 7f3a9e40dd7e79f842cc2c11fb8b548f9a802806, all input
sequences, generation manifests, native histories, truth, panel-specific seeds,
conditions, thresholds, clustering and phylogenetic settings. Preserve separate
fixed-length and variable-length analyses and their existing bootstrap seeds,
contrasts, multiplicity adjustment, failure reporting and no-zero-imputation
rules. Do not tune on corrected outcomes. Reusing exposed fixed-length seeds
does not turn that panel into independent confirmation.

The original source-only OrthoHMM checkout lacked compiled profile-alignment
code and silently built zero profiles. Retain those runs as defective-runtime
diagnostics, not estimates of the intended method. Correct inference uses the
separate CPU-native checkout `publication_method_native_v2` and frozen manifest
`publication_native_runtime_20260916.json`. This build also enables the compiled
search kernels; it is a documented runtime change, not merely a path rename.

## Execution And Reuse

- Freeze new manifests and use fresh output directories. Run only OrthoHMM
  high sensitivity and satellite_v2, with the original four-CPU allocation.
- Require the native manifest and actual profile smoke at preparation and
  execution. Each OrthoHMM process records verified runtime evidence before
  and after inference. Missing provenance stops execution before the tool runs.
- Keep the original OrthoFinder configurations and outputs. Never overwrite
  them or relaunch them as an implicit side effect of the OrthoHMM correction.
  Both successful and failed comparator outcomes remain eligible for audit;
  reuse itself does not certify completion or accuracy.
- Check dataset equality and identical scientific arguments while constructing
  reuse manifests. Freeze the original comparator manifest hash explicitly.
- Retain shared-machine timing caveats. Do not compare corrected-runtime
  OrthoHMM timings with defective-runtime timings as an optimization result.
- Require successful label-blind historical replay before submitting corrected
  simulation inference. The native smoke alone does not satisfy that gate.

## Scoring Gates

Both arrays must be terminal before mixed-provenance assembly. Supply and
verify each original array ID, executor revision, method manifest, input hashes,
output inventories and native completion evidence. Select only corrected
OrthoHMM rows and original OrthoFinder rows. Preserve their separate scheduler,
execution evidence, method-manifest hash, failures and resource records.

The current assembler rejects source-only OrthoHMM runs. Corrected OrthoHMM
requires per-process before/after native runtime checks and a passing profile
smoke in addition to the existing source, input and native-output gates.
Duplicate or mismatched datasets, differing truth, or unverified reused
execution provenance must not silently produce a comparative table.

After all gates pass, regenerate tables from corrected machine-readable
results and link the retained defective-runtime snapshots as audit evidence.
Do not overwrite those snapshots, pool runtime variants, or impute failures.

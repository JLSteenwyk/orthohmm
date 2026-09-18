# Corrected SwissTrees Primary Strata Execution

Status: implementation tested; real corrected outcomes not yet evaluated.
The unchanged [protocol](CORRECTED_SWISS_SEQUENCE_STRATA_PROTOCOL_20260918.md)
defines the estimand, resampling and multiplicity. This execution note does
not add endpoints or authorize tuning from outcomes.

## Evidence Inputs

`benchmark_tools/run_corrected_swiss_strata.py` requires:

- The complete eight-cell corrected factorial admission inventory and its
  SHA-256. It reruns the raw SwissTrees count audit before selecting cells
  p1c0r0 (high-sensitivity replay) and p1c1r1 (candidate expansion and native
  phylogenetic pairs). Other cells are not substituted for missing inputs.
- The full OrthoFinder3.1.5 corrected assessment admission and its SHA-256.
  It reruns the comparator raw-count audit; sequence-only output is rejected.
- The frozen baseline raw-count audit, used only to anchor reference
  identities and labels, never as a source of corrected prediction counts.
- The frozen `corrected_swiss_sequence_strata_20260918.json` and protocol.
  Input-only descriptors are reused and all bins are recomputed. The runner
  does not reread the staged FASTAs; descriptor provenance is the retained,
  hash-bound input-only extraction. Every selected count family must have
  exactly the frozen descriptor membership.

Seventeen analysis dependencies are hash-pinned. Source inputs and freshly
reconstructed evidence records are checked again after computation. Output
must be fresh. The result retains both reconstructed count reports, source
records, method bindings, missing bins and all27primary endpoints.

## Command Interface

Use the dedicated analysis environment with pinned NumPy and the frozen
executor revision containing this runner. The exact paths/hashes for the
first three arguments must come from terminal admitted results, not predicted
filenames or an unverified success tag:

```bash
python benchmark_tools/run_corrected_swiss_strata.py \
  --inventory "$FACTORIAL_ADMISSION_INVENTORY" \
  --inventory-sha256 "$FACTORIAL_INVENTORY_SHA256" \
  --orthofinder-admission "$FULL_OF_ASSESSMENT_ADMISSION" \
  --orthofinder-admission-sha256 "$FULL_OF_ADMISSION_SHA256" \
  --baseline benchmark_tools/results/qfo_swiss_counts_20260917.json \
  --strata benchmark_tools/results/corrected_swiss_sequence_strata_20260918.json \
  --protocol benchmark_tools/results/CORRECTED_SWISS_SEQUENCE_STRATA_PROTOCOL_20260918.md \
  --output "$FRESH_RESULT_PATH"
```

The runner reconstructs raw scorer evidence but does not rerun upstream
inference or scoring admissions. Its default100000draws and seed20260924
cannot be overridden from the CLI. Results remain development-exposed,
conditional family-bootstrap evidence, not independent confirmation.

## Validation And Remaining Work

Ninety-four focused tests passed, including18new driver tests and18kernel
tests. Driver tests use the real input-only stratum inventory but synthetic
prediction counts; orchestration tests mock the expensive raw auditors,
whose existing tests are included separately. They cover exact cell selection,
historical/mixed-method rejection, descriptor/count membership, recomputed
bins, changed sources, post-computation mutation and fresh-output guards.
The CLI help command also succeeds.

No corrected primary stratified result exists yet. Terminal corrected HMM,
phylogenetic and full OrthoFinder assessments remain prerequisites. The
all-method descriptive display, secondary overlapping strata, figure/table
integration and result-specific independent review remain to be completed.

# Saved FAS Complexity Exposure

This post hoc diagnostic intersects the saved corrected FAS samples with the
1,143 protein identifiers above the native path limit in the
[complete annotation panel](FAS_NATIVE_ANNOTATION_PATH_PANEL_20260928.md).
It preserves the separation between precomputed scores and newly calculated
scores and does not alter benchmark values.

## Results

The [generated table](FAS_SAVED_COMPLEXITY_TABLE_20260928.md) and
[full audit](fas_saved_complexity_exposure_20260928.json) show zero pairs
touching flagged proteins among saved newly calculated scores for every
method. Flagged proteins occur in precomputed saved pairs for all eight
methods: 52 for high-sensitivity OrthoHMM, 181 for satellite_v2, 30 for full
OrthoFinder, 14 for its sequence-only checkpoint, 27 for SonicParanoid,
1,633 for ProteinOrtho, 8 for FastOMA and 4,122 for recovered OrthoMCL.

The full lookup pass streamed 59,962,787 entries. It independently confirms
3,023,018 unique saved precomputed pairs and their values, and that 70,528
unique saved newly calculated pairs have no valid precomputed value. It
handles 2,762,101 relevant canonical overwrites from opposite accession orders.
These are unique-pair totals across methods, not independent observations.
All 46 focused tests pass. The run used ijson 3.5.0; parser-module hashes are
recorded along with the frozen lookup and all checked inputs.

## Checks

The diagnostic checks the retained annotation panel, its completeness review,
the corrected attrition audit and all saved raw score hashes. It rejects any
saved-score endpoint absent from the annotation inventory. A pair with two
flagged endpoints is counted once, and an empty category has no mean rather
than a zero score.

Initial stratum assignments use native raw-row order and logged stratum sizes.
An independent streaming join then checks those assignments against the frozen
precomputed JSON resource and reproduces every saved precomputed value within
1e-12. Opposite accession orders are canonicalized as in the native loader,
and the last valid canonical entry wins. Duplicate relevant JSON keys are
rejected because their insertion-order semantics would need separate handling.
An initial guard incorrectly treated distinct opposite-order keys as ambiguous;
the corrected checker and its tests explicitly model native overwrites.

The stream examines the full resource for membership of saved pairs; this is
not validation of every unused score, the underlying FAS algorithm or its
annotation correctness. Parser version and loaded parser-module hashes are
retained in the result. This is a read-only diagnostic with no rescoring.

## Interpretation

Precomputed values need not obey the current new-scoring complexity cutoff.
Their presence for flagged proteins is therefore not evidence of corrupt
scores and does not justify deleting them. Conversely, absence of flagged
proteins among saved new values would be consistent with the native omission
mechanism, not proof of which historical requested pairs were omitted.

The saved-sample intersection does not recount all eligible predicted pairs,
recover missing sample identities, quantify method-specific selection bias or
produce a family-aware comparison interval. Conditional category means are
descriptive diagnostics, not replacement FAS endpoints.

```sh
python -B -m benchmark_tools.audit_fas_complexity_exposure --output NEW_EXPOSURE.json
python -B -m pytest -q tests/unit/test_audit_fas_complexity_exposure.py \
  tests/unit/test_audit_fas_sample_attrition.py tests/unit/test_review_fas_path_panel.py
```

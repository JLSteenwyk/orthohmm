# All-Method Corrected GO/EC Pair Panel

The [machine-readable panel](qfo_scored_pair_panel_20260927.json) extends the
[three-method diagnostic](QFO_SCORED_PAIR_OVERLAP_20260927.md) to all eight
methods in the unchanged corrected comparison v7 manifest. All 56 pairwise
comparisons (28 per endpoint) have zero differing serialized scores on shared
pairs. Shared set sizes range from 57,090 to 230,397 pairs across comparisons.

| Method | GO scored pairs | EC scored pairs |
| --- | ---: | ---: |
| OrthoHMM high sensitivity | 145,619 | 185,664 |
| OrthoHMM phylogeny satellite_v2 | 84,211 | 117,460 |
| OrthoFinder 3.1.5 full | 163,557 | 175,361 |
| OrthoFinder 3.1.5 sequence-only checkpoint | 1,979,084 | 1,660,940 |
| SonicParanoid 2.0.9 | 197,930 | 244,356 |
| ProteinOrtho 6.3.6 | 89,609 | 143,135 |
| FastOMA 0.3.5 | 171,132 | 178,598 |
| OrthoMCL 1.4 | 175,181 | 220,622 |

All 16 raw counts match admitted assessed-relation counts, and all rounded
means match admitted means within 5.051e-7. Forty-three source/input record
entries were verified before/after analysis. The full manifest hash fixes
method identity and order; non-admitted, duplicate or changed method panels
are rejected. FastOMA's supplied-tree and OrthoFinder checkpoint scopes remain
unchanged. These counts concern scored eligible relations, not all predictions.

At retained six-decimal precision, differences among these methods' GO/EC
means reflect different pair sets and denominators. Intersection-only means
would erase this distinction and are not substitute endpoints. Identical
rounded pair scores do not prove identical full-precision scores, annotation
provenance or scorer internals. This arithmetic result supplies no independent
resampling units, confidence intervals, causal explanation or accuracy ranking.

## Reproduction

```bash
python -m benchmark_tools.run_qfo_scored_pair_panel --repo . \
  --output NEW_PANEL.json
python -m pytest -q tests/unit/test_run_qfo_scored_pair_panel.py \
  tests/unit/test_compare_qfo_scored_pairs.py tests/unit/test_audit_qfo_go_ec.py
```

The panel writes only a fresh final report after all checks pass. No native
inference/scoring is rerun. All 35 focused tests pass, including complete-panel
enumeration, changed manifest, duplicate/non-admitted methods, count/mean
disagreement and overwrite refusal. Retained raw inputs are needed to execute
the batch; this report is not a portable raw-data release or publication gate
closure.

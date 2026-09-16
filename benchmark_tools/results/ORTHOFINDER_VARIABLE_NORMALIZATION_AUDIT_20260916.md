# Variable-Length OrthoFinder Normalization Audit

## Scope

The [machine-readable diagnostic](orthofinder_variable_normalization_diagnostic_20260916.json)
recomputes all 640 ordered species-pair score matrices across all ten divergent
variable-length seeds, including five valid and five failed runs. It uses the
frozen OrthoFinder 3.1.5 interpreter, package inventory, native BLAST parser and
normalization functions. Existing output inventories are checked before and
after. No installed code, input, prediction, score, or exclusion is changed.

## Reproduced Failure

Nonfinite normalized scores occur in exactly the five seeds with nonfinite
saved OrthoFinder graphs. None occurs in the five seeds whose graphs passed.
All seven affected matrices are within-species searches, with two non-self
stored hits at a single length product. Their two-parameter log-linear fits
have design rank one. The fitted intercepts range from approximately -342 to
-404; the native implementation separately exponentiates that intercept,
overflowing its scale factor and producing nonfinite normalized values.

| Seed | Native species pair | Length product | Raw non-self hits | Nonfinite normalized entries |
| --- | --- | ---: | ---: | ---: |
| 20261101 | 1,1 | 151321 | 2 | 2 |
| 20261101 | 5,5 | 151321 | 2 | 2 |
| 20261102 | 2,2 | 53824 | 2 | 2 |
| 20261107 | 7,7 | 43681 | 2 | 2 |
| 20261108 | 1,1 | 148225 | 2 | 2 |
| 20261109 | 4,4 | 193600 | 2 | 2 |
| 20261109 | 6,6 | 233289 | 2 | 2 |

Species indices are OrthoFinder's internal indices. Identical-sequence
self-hits are excluded by the native parser; these are not missing-self-hit
errors. All proteomes have 82-92 distinct protein lengths. Heterogeneous
proteome-wide lengths therefore do not guarantee sufficient variation in
each local fit.

Across the 640 matrices, 23 fitted subsets have rank one with at least two
observations: seven yield nonfinite values, three yield no stored normalized
scores, and thirteen yield finite stored values. The three empty outcomes
occur in already-failed seeds 20261108/20261109, with large positive fitted
intercepts. Rank deficiency alone is not identical to the observed nonfinite
failure condition. No normalization call raises an exception in this audit.

## Interpretation

This establishes a reproducible local numerical failure in the frozen
implementation, distinct from the original all-proteins-equal-length panel.
It does not establish failure frequency on real proteomes or justify changing
the comparator retrospectively. Keep the five existing exclusions and
complete-case confidence intervals as reported. Finite native outputs are not
newly excluded on the basis of this post-outcome diagnostic. Any alternative
normalizer would be a separately specified comparator configuration, not an
undocumented repair to the published baseline. No upstream issue was filed.

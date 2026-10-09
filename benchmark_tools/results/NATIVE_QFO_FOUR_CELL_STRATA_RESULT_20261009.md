# Four-Cell Native SwissTrees Fixed-Stratum Results

The separately named prospective correction was committed/pushed at
`67dc83ab` after 137 focused/adjacent tests passed. The original v1 attempt
[failed before writing output](native_qfo_four_cell_strata_failure_20261009_v1.json)
because three retained distance rows used `native_pair`, whereas the consumer
expected `resolved_native_pairs`. The [compatibility amendment](NATIVE_QFO_FOUR_CELL_COMPATIBILITY_AMENDMENT_20261009.md)
authorizes a copied metadata view only, preserving original labels and pins.
Old sources, inputs, tests and the scientific protocol remain unchanged.

## Verified Output

The [v2 report](native_qfo_four_cell_strata_20261009_v2/report.json),
125,996 bytes, SHA256
`ef3b9de78197a946399543b5aec0b35f0cbe9c89cce643da1a7024cae48fda73`,
contains 72 cell/family rows, 92 cell/bin scores, 69 contrast/bin differences,
and all three original-to-normalized metadata mappings. All 69 old bin scores
and 46 old differences reproduce within absolute 1e-12; inherited distance
differences are explicitly converted from percentage points to raw units.

The independently executed [rational readback](native_qfo_four_cell_strata_readback_20261009_v2.json),
8,804 bytes, SHA256
`1260c1c87ec4d3cbaa606a52225e94fc1fa084a689ba8af9e721f2f14966d1c2`,
checks the truthful v2 schema and original metadata mapping separately before
reusing the unchanged independent Fraction/table kernel. It checks all counts,
statistics, 23 human table rows, TSV values and 26 direct references. Both
actual commands exited zero; the [execution record](native_qfo_four_cell_strata_execution_20261009_v2.json)
retains their code, commands and tool output. This is retained-record validation,
not repeated raw scoring or transitive scientific admission.

## Profile Contrast

The [complete table](native_qfo_four_cell_strata_20261009_v2/TABLE.md) reports
all 23 fixed bins and all three conditional contrasts, including five empty
bins as NA. Scores are percent and differences are percentage points; JSON
and TSV values retain raw units. Initial HMM search is on in every cell.
P1/C0/R1 minus P0/C0/R1 measures downstream profile refinement at C0/R1,
not total HMM benefit or a profile effect at other settings.

Overall macro F1 changes by -0.320215 percentage points. Lower-entropy families
have +0.132851 points, while higher-entropy families have -0.762755. Families
above the fixed median model distance have +0.160465 points, while those at
or below it have -0.794278. Each of these bins contains nine families. These
contrasts separate CASP and GH14, not many independent profile successes:
CASP loses three false positives (FP to TN), and GH14 loses one true positive
(TP to FN); the other 16 families' counts are unchanged.

Both changed families are in the not-short-relative bin, lower median Pfam
type-count bin, lower repeat-fraction bin and upper reference-duplication bin;
their F1 contrasts are -0.499153, -0.495776, -0.377910 and -0.679748 points,
respectively. The complementary bins have exactly unchanged counts and scores.
All-family, no-explicit-fragment and not-concentrated rows repeat the overall
contrast and are not independent findings. Empty bins do not support any effect.

The existing [complete pair localization](NATIVE_PROFILE_LOCALIZATION_RESULT_20261007.md)
places all four removed calls at pre-reconciliation candidate separation.
These retained counts and memberships explain the descriptive bin pattern;
they do not identify a causal HMM edge, validate the gene/species trees or
demonstrate a general dependence of performance on entropy or divergence.

## Limits

No subgroup confidence intervals, new bootstrap, cutoff selection, method
change or default promotion occurred. Existing aggregate intervals are not
subgroup intervals. The fixed bins are retrospective, overlapping and
development-exposed. Sequence flags are not literal fragment truth, Pfam
descriptors are incomplete, and reference duplication/model distances are
not known ancestral histories or dated divergence. Missing fresh factorial
cells remain missing. Recovered timing remains failed/ineligible; shared-host
timing distortion is unknown and potentially tool-dependent. General
superiority over full OrthoFinder and publication readiness are not established.

The remaining integration task is to add this complete profile-bin extension
to a new machine-generated manuscript without changing frozen parent drafts
or relabeling the older three-cell figure as a four-cell result. Other original
scientific limitations remain open; no new archive candidate is warranted.

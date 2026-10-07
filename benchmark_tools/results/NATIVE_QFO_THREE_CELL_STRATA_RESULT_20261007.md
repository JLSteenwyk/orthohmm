# Native SwissTrees Three-Cell Fixed Strata

The [complete table](native_qfo_three_cell_strata_20261007_v1/TABLE.md),
[report](native_qfo_three_cell_strata_20261007_v1/report.json) and
[independent rational readback](native_qfo_three_cell_strata_readback_20261007_v2.json)
extend the unchanged two-cell sequence/domain/duplication bins to admitted
candidate-expansion cell8. All18families/563canonical reference accessions
match; all40prior cell/bin scores reproduce, all54family rows/60score rows/
40conditional contrast rows are checked. Empty bins remain NULL/NA. The
repeated all-family rows and overlapping bins are not independent findings.

The [protocol](NATIVE_QFO_THREE_CELL_STRATA_PROTOCOL_20261007.md) was pushed
762cbdb0 before selected subgroup calculation; tested sources a5952bbc were
pushed before export. The original reader encountered a string-path bug,
preserved in its [failure receipt](native_qfo_three_cell_strata_readback_failure_20261007_v1.json).
The [versioned amendment](NATIVE_QFO_THREE_CELL_STRATA_READER_AMENDMENT_20261007.md),
new reader/tests and failure receipt were pushed d8b0981e before a fresh
v2readback. The export/report were not rewritten or rerun. Both actual
successful postprocessing runs exit0 in isolated scientificPython3.10.13;
see [execution provenance](native_qfo_three_cell_strata_execution_20261007_v1.json).
Their0.05s durations/RSS are not native inference costs. No job or scientific
score/timing admission changed; cell7's failed timing remains ineligible.

## Interpretation

All three cells have initial HMM search on and downstream profile refinement
off. R1 and C1 are separate conditional contrasts against the same baseline,
not candidate/reconciliation interaction or head-to-head OrthoFinder tests.
The statistic is macro family precision/recall then harmonic F1, using
the original halved integer counts plus1 prior, not pooled pairs or
mean-family F1. Changes below are percentage points, not relative percent.

Across all18families, C1 decreases precision3.523pp and increases recall4.375pp,
with F1 changing -0.347pp. Its conditional F1 changes are:

| Frozen bin pair | Families | First-bin C F1 change | Second-bin C F1 change |
|---|---:|---:|---:|
| Higher/lower entropy | 9/9 | -0.736 | +0.599 |
| Short-relative/not-short-relative | 7/11 | -0.703 | -0.236 |
| Median Pfam types at least2/below2 | 6/12 | +0.643 | -0.672 |
| Repeated-Pfam fraction at least1/4/below1/4 | 3/15 | +0.987 | -0.591 |
| Lower/upper reference duplication fraction | 9/9 | -0.467 | -0.273 |

Thus the aggregate trade-off is not uniform: the lower-entropy, higher-Pfam
and repeated-Pfam bins show small descriptive F1 increases, while the other
listed bins decrease. Candidate precision decreases and recall increases in
every nonempty bin. The three-family repeated-Pfam bin is particularly small;
these descriptions do not establish a reliable rule for choosing merges.
The full table also preserves R1's prior precision/recall/F1 results. No
subgroup cutoff or default is selected, and no new bootstrap/interval is
attached to these bins. Retained aggregate intervals belong to their original
aggregate contrast, not subgroup inference; neither aggregate F1 contrast
establishes superiority from its adjusted interval.

Length/entropy are proxies, not calibrated divergence or fragment truth.
Pfam descriptors are annotations, not complete architecture validation.
Reference informative-node duplication fractions are not ancestral histories;
default S was not explicit speciation. These development-exposed associations
do not identify causal biological mechanisms or independent generalization.
Source/raw feature extraction and prior readers were not repeated: the new
reader independently checks bounded artifacts and exact arithmetic, not
full transitive scientific admission. Existing main text and rc5 package
remain unchanged; this is a new error-analysis supplement, not a regenerated
archive or a declaration that the full publication goal is complete.

# Fixed-Candidate Gene-Tree Oracle Results

## Methods

This development-exposed mechanism diagnostic retains all70 generating-
species-tree simulation cells, seven conditions and ten seeds each. It
holds candidates, species trees, frozen reconciliation rules and constraint
policy fixed. Only actually inferred single-ancestral-family candidates
receive generating gene trees, either with their generating root or with
frozen minimum-duplication/loss rerooting. All bypasses/ineligible candidates
retain native predictions. The [prospective protocol](SIMULATION_GENE_TREE_ORACLE_PROTOCOL_20261004.md)
and tested runner were committed/pushed before new scoring.

## Accuracy

F1 values are percentages, arithmetic means of ten dataset scores; changes
are percentage points. These are descriptive finite-panel effects, not
population confidence intervals, significance tests or independent validation.

| Condition | Inferred | Generating root | Generating, rerooted | Generating root - inferred |
| --- | ---: | ---: | ---: | ---: |
| baseline | 99.493 | 99.699 | 99.679 | +0.206 |
| divergent | 69.631 | 69.768 | 69.731 | +0.136 |
| divergent_turnover | 71.672 | 72.168 | 72.029 | +0.496 |
| missing20 | 99.403 | 99.622 | 99.578 | +0.219 |
| taxon_count_control | 99.316 | 99.554 | 99.536 | +0.238 |
| turnover | 98.990 | 99.724 | 99.441 | +0.734 |
| uneven_taxa | 99.028 | 99.203 | 99.203 | +0.174 |

All70 baseline pair sets reproduce native output; independent readback also
matches every original retained generating-species-tree score. There are
10,125 candidates and 1,377 eligible gene-tree controls.
Among eligible controls, 300 have unrooted topology disagreement;
30 more differ only in root (330 rooted disagreements total).
This does not establish that every disagreement causes an accuracy error.

## Residual Errors

Counts below sum ten seeds per condition; they are not the denominator of
the macro-mean F1 above. The cross-candidate count is independently derived
from truth and the native candidate partition, not inferred from a score.

| Condition | Generating-root FN | True pairs across candidates | Across-candidate share | Eligible within-candidate FN |
| --- | ---: | ---: | ---: | ---: |
| baseline | 177 | 177 | 100.000% | 0 |
| divergent | 12,766 | 12,736 | 99.765% | 30 |
| divergent_turnover | 13,104 | 13,072 | 99.756% | 32 |
| missing20 | 138 | 138 | 100.000% | 0 |
| taxon_count_control | 92 | 92 | 100.000% | 0 |
| turnover | 196 | 196 | 100.000% | 0 |
| uneven_taxa | 162 | 162 | 100.000% | 0 |

In both divergent conditions more than99.7% of residual generating-root
false negatives are fixed upstream by different candidate membership.
Replacing gene trees alone cannot recover those pairs. This localizes this
deficit upstream of gene-tree reconciliation, but does not distinguish HMM
search, grouping and candidate expansion. Generating trees are not uniformly
perfect: retain the within-candidate errors, bypass false positives and
losses as well as gains. No defaults or benchmark scores are changed.

## Limits And Reproduction

True gene trees/roots are unavailable in ordinary inference. The oracle
also changes lengths/support and is not topology-only causal evidence or an
accuracy upper bound. Native pair scoring does not validate ancestral-copy
root-HOG membership. This narrow mechanism result does not close all error-
analysis, independent-generalization or uncertainty requirements.

The [compact readback](simulation_gene_tree_oracle_readback_20261004.json)
retains all70 score/decomposition rows; the detailed local report records
every candidate and all source/input identities. Readback independently
checks counts/partitions and original scores, not the new oracle pair sets.
It rehashes 6,334 input/source identities.

Detailed local report: `/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmark_tools/results/simulation_gene_tree_oracle_20261004.json`.
SHA256: `d2299d73238f1a4bb511320add9a62dbb427bc1ea9b4f414e74aefab20e7ad53`.

The [execution receipt](simulation_gene_tree_oracle_execution_20261004.json)
records successful execution, tests and the failed first readback. The first
readback9e450025 rejected before output because it expected `scored` rather
than the retained `complete` status. Corrected9893f6a0 explicitly validates
all70 baseline completions; no seed or outcome was dropped. The44 diagnostic/
readback tests and two renderer tests pass. No completed inference or timing
measurement was repeated.

This is a small read-only diagnostic on the shared Threadripper, not a
matched timing run. No unrelated analysis was stopped or modified. Existing
rc2/PDF artifacts remain unchanged; this is a later scientific addendum.

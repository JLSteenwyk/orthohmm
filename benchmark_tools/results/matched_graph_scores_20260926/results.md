# Matched-Recall Graph Control

All values below are percentages or percentage-point differences.
Intervals are the prespecified Bonferroni-adjusted eight-contrast F1 intervals.

| Condition | HMM F1 | DIAMOND F1 | Difference (pp) | Adjusted interval (pp) | Wins/ties/losses |
|---|---:|---:|---:|---|---|
| baseline | 95.7240 | 95.3749 | +0.3492 | [+0.0000, +0.8995] | 2/3/0 |
| divergent | 59.0789 | 49.4565 | +9.6223 | [+6.0348, +13.7165] | 5/0/0 |
| divergent_turnover | 59.5801 | 48.8969 | +10.6832 | [+5.0577, +17.5016] | 5/0/0 |
| missing20 | 95.6338 | 95.3321 | +0.3017 | [+0.0000, +0.6385] | 3/2/0 |
| taxon_count_control | 95.5191 | 95.1637 | +0.3553 | [+0.0000, +1.1295] | 3/2/0 |
| turnover | 84.9816 | 84.5499 | +0.4317 | [+0.0000, +1.4301] | 3/2/0 |
| uneven_taxa | 94.7504 | 94.8518 | -0.1014 | [-1.5113, +0.9457] | 2/2/1 |
| overall | 83.6097 | 80.5180 | +3.0917 | [+1.6782, +4.5313] | 5/0/0 |

## Limits

- Development-exposed simulations; not independent confirmation or real-data sensitivity matching.
- Five paired seed blocks give limited tail resolution; percentile intervals have approximate coverage.
- Initial-search graph control only: profile expansion, candidate expansion and phylogeny are off.
- Cluster-derived cross-species pairs are not native reconciled ortholog predictions.
- Matching search recall does not equalize score distributions, computational effort or hit identities.
- Shared-host incremental resources are descriptive, not controlled efficiency evidence.

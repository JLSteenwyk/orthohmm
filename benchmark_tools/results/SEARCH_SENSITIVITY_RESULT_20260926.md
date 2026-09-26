# Simulation Search-Sensitivity Result

All 70 cells in Slurm array 22180 completed with exit 0. The independent reader
validated frozen input/truth identities, installed runtime/source checks,
commands, environment, native output inventories, HMM species-pair counts,
gene identities, DIAMOND lengths, finite scores and duplicate-free directed
hits. No inference retries were used. Raw output remains under
`benchmarks/work/search_sensitivity_panel_20260926`; every execution receipt
and its output hashes are retained.

[Machine-readable results](search_sensitivity_results_20260926.json), 349,785
bytes, SHA-256
`d58b1603e563eda27d90cc113950fb28091c94d57233a60ec27b32e297ad729c`.
[Frozen protocol](SEARCH_SENSITIVITY_CALIBRATION_PROTOCOL_20260926.md).
The scorer was committed as `e7076e1` before execution; 64 focused tests pass.
A separate NumPy aggregation of reported per-dataset recalls agrees with the
exact-fraction cutoff selection and overall reporting difference within 1e-12.

## Calibration

The 35 calibration datasets have fixed HMM primary recall **83.8449%**.
The prespecified minimum absolute recall-gap rule selects DIAMOND **E <= 1e-40**
(83.3844%). The complete calibration grid is retained below and in the JSON.
The HMM search retains its native strict E < 1e-4 rule.

| DIAMOND cutoff | Calibration recall (%) |
|---|---:|
| 1e-100 | 56.3769 |
| 1e-80 | 65.9647 |
| 1e-60 | 74.7861 |
| **1e-40** | **83.3844** |
| 1e-30 | 86.8808 |
| 1e-20 | 89.7165 |
| 1e-15 | 90.9118 |
| 1e-10 | 91.8695 |
| 1e-8 | 92.1947 |
| 1e-6 | 92.4832 |
| 1e-4 | 92.6759 |
| 1e-3 | 92.7090 |
| 1e-2 | 92.7216 |
| 1e-1 | 92.7233 |
| 1 | 92.7249 |

## Reporting Split

Recall gives equal weight to eligible ancestral families within each dataset,
then equal weight to datasets. Each condition has five reporting seeds.
Differences are DIAMOND minus HMM in percentage points.

| Condition | HMM recall (%) | DIAMOND recall (%) | Difference (pp) |
|---|---:|---:|---:|
| Baseline | 95.7334 | 96.5119 | +0.7786 |
| Divergent | 30.9157 | 30.4699 | -0.4459 |
| Divergent + turnover | 30.5774 | 29.8276 | -0.7497 |
| Missing 20% | 95.8910 | 96.3502 | +0.4592 |
| Taxon-count control | 95.5266 | 96.0996 | +0.5730 |
| Turnover | 94.7206 | 94.9050 | +0.1844 |
| Uneven taxa | 94.2584 | 95.4688 | +1.2103 |
| **Overall** | **76.8033** | **77.0904** | **+0.2871** |

The frozen gate **passes**: overall absolute difference is below 2 percentage
points and each condition is below 5 points. Reporting datasets contain 204,954
directed cross-species ancestral-homology pairs. HMM retains 135,644
cross-species hits (135,642 homologs and 2 nonhomologs); selected DIAMOND retains
135,706 (all homologs). These pooled counts are not the primary mean statistic
and are not independent observations.

## Interpretation and Next Step

This establishes a bounded matched-recall search control on this particular
simulation panel. It does not establish real-proteome sensitivity equivalence,
equal computational effort, orthology F1 improvement or general superiority.
Ancestral homology includes paralogs. These simulations were previously exposed
to development; the seed split is not independent validation of the method.
Conditions sharing seeds share histories. No confidence interval or significance
claim is made. DIAMOND's broad native search was post-filtered using printed
E-values, not rerun at the selected threshold; shared-host timings remain
descriptive only.

A downstream experiment must now be separately prespecified: reuse these
frozen search outputs with identical graph construction/clustering and evaluate
orthology outcomes without further threshold selection. Do not transfer the
selected cutoff to QfO or OrthoBench as a claim of matched real-data sensitivity.
Existing real-data controls, scientific defaults and publication claims remain
unchanged. Dedicated timing and the other publication requirements remain open.

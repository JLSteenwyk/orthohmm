# Corrected VGNC Single-Block Sensitivity

Protocol and runner were committed and pushed as `3d49a44c` before calculation.
All eight methods and all 16,844 reference blocks are retained, producing
134,752 individual deletions and seven matched method contrasts against full
OrthoFinder. A deletion removes every incident scored cell once, including
cross-block false positives, then recomputes ratios from remaining counts.

| Candidate minus full OrthoFinder | Original F1 difference (pp) | Range after any one block deletion (pp) |
| --- | ---: | --- |
| OrthoHMM high sensitivity | -32.2613 | [-32.3167, -32.2344] |
| OrthoHMM phylogenetic | -8.6856 | [-8.7353, -8.6729] |
| OrthoFinder sequence-only | -84.5791 | [-84.6343, -84.5590] |
| SonicParanoid | -0.5752 | [-0.5877, -0.5624] |
| ProteinOrtho | -3.3650 | [-3.4193, -3.2737] |
| FastOMA | -3.8450 | [-3.8704, -3.8260] |
| OrthoMCL | -34.6520 | [-34.6910, -34.4040] |

Every contrast remains negative for all 16,844 deletions. In this narrow
fixed-table diagnostic, no single reference block accounts for a ranking
reversal. For phylogenetic OrthoHMM the minimum difference occurs after deleting
KRT31 and the maximum after deleting ABCG2. These extrema are descriptive
selections, not prespecified biological case studies or causal explanations.

**These ranges are not confidence intervals.** The deletions share scored
pairs, change the table and do not regenerate native eligibility or predictions.
Stability against one deletion does not establish stability against joint
deletions, new datasets or population sampling. The failed rare-error/shared-
clade uncertainty diagnostics remain unresolved. No benchmark point estimate,
method configuration or native statistical claim is replaced.

Seventeen focused influence/component tests pass. An independent standard-
library accumulation from each sparse cell table, followed by exact rational
TP/FP/FN ratios, checked all 134,752 exported deletion rows and all seven paired
contrast ranges and sign counts. Maximum discrepancy was 5.56e-17 in 0-to-1
units. Source and input records were checked before and after both calculations.

- [Full result and input records](corrected_vgnc_influence_20260928.json)
- [Independent arithmetic check](corrected_vgnc_influence_crosscheck_20260928.json)
- [Frozen exploratory protocol](CORRECTED_VGNC_INFLUENCE_PROTOCOL_20260928.md)

All deletion rows are retained locally at
`benchmarks/work/corrected_vgnc_influence_20260928/all_block_deletions.tsv.gz`;
their checksum is in the result. This is not a standalone archival bundle.

```sh
python -B -m benchmark_tools.diagnose_corrected_vgnc_influence \
  --source benchmark_tools/results/corrected_vgnc_blocks_20260926.json \
  --output /tmp/corrected-vgnc-influence-new
```

Use a fresh output directory and the retained local inputs. No native inference
or scoring run was repeated. Controlled timing and publication gates stay open.

## Executable Independent Check

The initial independent calculation is now available as a standard-library
command, separate from the primary incident-counting and scoring helpers:

```sh
python3 -B -m benchmark_tools.check_corrected_vgnc_influence \
  --mapping benchmark_tools/results/corrected_vgnc_blocks_20260926.json \
  --report benchmark_tools/results/corrected_vgnc_influence_20260928.json \
  --output /tmp/corrected-vgnc-rational-check-new.json
```

The [system-Python replay](corrected_vgnc_executable_replay_20260929.json)
checks all 134,752 deletion rows, full method counts/ratios, all seven paired
ranges/sign counts and the values at reported extremal blocks. It additionally
checks each saved F1 change against exact rational subtraction; the expanded
check has maximum error 1.12e-16 in 0-to-1 units. Twenty focused checker/runner
tests pass, including nonfinite values and malformed counts/cells. It does not
validate the top-ten influence ordering, stochastic assumptions, raw prediction
generation or independent cross-host restoration. Local sparse tables remain
required. Original results and the first independent receipt are unchanged.

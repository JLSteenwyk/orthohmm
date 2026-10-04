# Native Costs Linked To OrthoBench Configurations

All repeats retained; shared-host observations with unknown, potentially method-dependent contention.
Different private runtime from original factorial; native inference/output writing only, not whole workflow.

| Cell | Repeats | Exact partitions | Wall median [range], s | CPU median [range], s | Lifetime peak median [range], GiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| p1_c0_r0 | 3 | 3/3 | 2720.077 [2621.498, 3027.061] | 78335.831 [75126.859, 86415.319] | 10.555 [10.533, 10.730] |
| p1_c1_r1 | 3 | 1/3 | 4392.872 [4042.426, 4809.292] | 125264.737 [116000.996, 136293.332] | 10.511 [10.503, 10.603] |

## Partition Readback

| Run | Cell | Partition equal | Genes in changed groups | Reference genes in changed groups |
| ---: | --- | --- | ---: | ---: |
| 7 | p1_c0_r0 | True | 0 | 0 |
| 8 | p1_c1_r1 | True | 0 | 0 |
| 15 | p1_c0_r0 | True | 0 | 0 |
| 16 | p1_c1_r1 | False | 149 | 0 |
| 24 | p1_c1_r1 | False | 172 | 0 |
| 26 | p1_c0_r0 | True | 0 | 0 |

## Limits

- All six points are retained shared-host native CLI observations. Contention distortion is unknown and potentially method dependent; do not infer isolated tool speed or causal HMM/phylogeny overhead.
- Native launch-to-exit includes inference and output/metrics writing; excludes preparation, harness hashing, conversion and scoring. CPU includes wrapper work; lifetime cgroup memory peak includes launcher work.
- Frozen scientific source and original twelve-proteome bytes are checked, but the patched private runtime differs from the original factorial deployment. These are not costs of the original cached-stage executions.
- Phylogenetic candidate partitions differ despite identical prescribed settings and counts; two final partitions differ in nonreference groups. Retain every repeat, not only matching outputs. These costs describe the prescribed native configurations, not exact reproduction of every original factorial intermediate. Reference-touching final group equality is a structural check, not newly recomputed accuracy or proof of exact whole-output reproducibility.
- Native pair files, inferred trees and alignment bytes are not compared. No pairwise-orthology equivalence or cause of output differences is established.
- Reuse pinned terminal environment/resource/runtime reviews; do not repeat their raw forensic audits. New candidate files are checked as current readbacks, not independently attested historical write events.
- Fourteen other factorial configurations have no associated full-native costs in this report. Do not impute them, combine old/new memory scopes, or pool observations across deployments.

# Controlled Fragment Inputs Prepared

The prospective source milestones were committed/pushed before production use:
c5c81c1f (transform, scoring and retained controls),3456a889 (explicit runtime
inventory amendment),67b6b34d (sequential runner and actual prospective runtime).
Final affected unit suite:167 passed in2.39s; shell syntax check passed.

The single preparation executed successfully, with10 datasets and30 genuinely
new native inference identities. Prepared manifest:
`benchmarks/work/controlled_fragment_observations_20261009_v1/manifest.json`,
600745 bytes, SHA256
`8e58830911f60de5aec9ad616e5c02c12229702780e26216aa1371017637666d`.
Only this manifest is retained in Git; generated FASTA files and per-gene
coordinate records stay at its pinned local paths, not committed as raw data.

The baseline check revalidated all40 native admissions and original scores.
The first strict historical-runtime check failed before inference because the
OrthoHMM package inventory changed. The separate prospective runtime passed,
retaining all72 inventory differences and unchanged scientific dependencies:
`controlled_fragment_execution_runtime_20261009_v1.json`, SHA256
`4ea5d10bccd02ef42e49194908212d07af701c109db878e176f56c298c2ca4da`.
No historical inventory or output-equivalence claim is made.

The new runner's check-only execution passed against both actual hashes.
A separate inline reader imported none of the new transform, command or
validator modules. Using Bio.SeqIO and direct SHA256 ranking/substrings, it
checked every written observation, coordinate field, original/observed hash,
gene ID and species owner, exact parent-truth bytes, all30 command path changes
and unchanged scientific arguments, and the diagnostic checkpoint binding.

| Seed | Retained Genes | Fragment Observations |
| --- | ---: | ---: |
| 20261101 | 806 | 161 |
| 20261102 | 822 | 164 |
| 20261103 | 813 | 162 |
| 20261104 | 811 | 162 |
| 20261105 | 806 | 161 |
| 20261106 | 847 | 169 |
| 20261107 | 804 | 160 |
| 20261108 | 816 | 163 |
| 20261109 | 852 | 170 |
| 20261110 | 823 | 164 |
| Total | 8200 | 1636 |

All other sequences are unchanged. This is the prespecified synthetic
observation test, not natural fragment truth or independent biological
validation. No fragment inference or accuracy result existed at preparation.
Use the retained protocol for scoring and failures; do not tune on new scores.

Execution uses the existing Threadripper scheduler,16CPUs/16GiB, ascending
seeds and high-sensitivity/satellite_v2/full-OrthoFinder order within each seed.
Shared-host timings remain descriptive; cached baseline costs are not paired
runtime controls. No quiet window, DGX, automatic retry or baseline rerun.

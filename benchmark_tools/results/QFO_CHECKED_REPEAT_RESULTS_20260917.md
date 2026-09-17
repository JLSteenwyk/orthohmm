# Checked QfO Initial-Graph Repeats

Slurm job **21328** completed with exit **0:0** in **38:14**. Executor
`f6ad87c` ran exactly three fresh, sequential, single-CPU workers. Each used
streamed Python-integer pairs to construct the same saved graph, followed by the
unchanged frozen CPM optimizer (resolution0.1, seed4, two iterations).

## Results

| Repeat | Worker wall seconds | Groups | Complete unique gene coverage | Same partition as repeat0 |
| --- | ---: | ---: | ---: | --- |
| 0 | 1044.452333 | 349898 | 976504 | Yes, self-check |
| 1 | 751.549762 | 349898 | 976504 | Yes, byte-identical |
| 2 | 463.853209 | 349898 | 976504 | Yes, byte-identical |

All three pairwise comparisons are byte-identical and membership-identical.
Every worker recorded the complete saved endpoint and weight fingerprints both
before and after optimization, with24,148,515 edges and976,504 vertices. The
constructor input's oriented int32 bytes also match a fresh saved-array hash.
No failed run was retried, and no accuracy labels or scores selected a partition.

Common partition SHA256:
`8c162782bba955078a8f7e4e6088df5b34c9f15b055803224e513fdea8b33cdd`.
Common canonical ordered endpoint SHA256:
`dd9c0c07cc9d9891d0ce3b50425a9b9c867caf0fe4d0400025425e0b925d4cd5`.
Common ordered weight SHA256:
`d398011be5171a8b6d8f41eab6972a081de5e9419f69536ec69494c746b4ecbf`.

## Independent Admission

`admit_qfo_checked_repeats.py` verified completed scheduler accounting, the exact
three-repeat inventory and frozen executor, commands, native runtime, preserved
worker files and parent-report agreement. It independently reconstructed saved
graph hashes, validated optimizer/constructor observations and full unique gene
coverage, recomputed all pairwise comparisons, and checked308 provenance records.

Raw report: `benchmarks/results/qfo_checked_repeats_v1/results.json`, SHA256
`26e3bee1cc869be5df02f6160aba6fb2bb024b98a98eb348b8a5d29af12dbad5`.
Admitted snapshot: [qfo_checked_repeats_verified_20260917.json](qfo_checked_repeats_verified_20260917.json), SHA256
`78f51f5ce703caf39307c5e518ad737acdf655dbe0926a34f2c20bbdec6d03f1`.

## Interpretation and Next Gate

This is bounded evidence for intact graph construction and matching optimization
outputs under this checked execution path. It does not prove general determinism,
a specific native-library defect, or equivalence to historical complete outputs.
The audit validates preserved observations, not retrospective live-memory access.
Input-format conversion and instrumentation change allocation and timing.

Worker times are shared-machine diagnostic costs, including validation overhead;
their decline is not a measured speedup. Scheduler MaxRSS is absent and is not
filled with an inferred value. No benchmark accuracy estimate has changed.

Next, preserve and check every clustering payload in one complete cached QfO
replay, including singleton/multipass and profile stages, without selecting among
retries. Require actual nonzero profile construction where expected, complete
output coverage, and input/runtime identity before admitting ablation inputs.
Report historical-stage comparisons without requiring agreement as a condition
for retaining an otherwise valid run. Any difference remains explicit and needs
investigation; do not silently substitute the new output for a historical baseline.

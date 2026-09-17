# Candidate Neighborhood Preparation Admission

Completed preparation job21314 independently verified, without accuracy scoring.
The unchanged control reproduces both frozen candidate-partition and merge-trace
bytes. All five arms retain all251,378 input proteins from62,885 HMM seed groups.

| Arm | Minimum normalized support | Minimum margin | Candidate groups | Reconstructed merges |
| --- | ---: | ---: | ---: | ---: |
| Control | 0.030 | 1.5 | 54,445 | 8,440 |
| norm_low | 0.024 | 1.5 | 54,370 | 8,515 |
| norm_high | 0.036 | 1.5 | 54,540 | 8,345 |
| margin_low | 0.030 | 1.2 | 52,912 | 9,973 |
| margin_high | 0.030 | 1.8 | 55,434 | 7,451 |

Validator checks terminal scheduler status, frozen manifest hashes, exact arm
order, applied overrides versus the separately retained nominal wrapper report,
one engine invocation per arm, all70 provenance records before and after,
complete/disjoint input coverage, constraint membership, and reconstruction of
each partition from its seed groups and ordered merge events. Count summaries
are independently recomputed from native files.

This admission does not independently recompute search evidence or verify that
each threshold decision was biologically correct. It establishes integrity for
downstream evaluation, not accuracy or general robustness. No variant was chosen
using reference scores; none is promoted as a default.

Next: run the four changed candidate partitions through the frozen inferred-tree
pipeline with their own constraints. Existing raw gene-tree checkpoints may be
reused only when the native checkpoint validation establishes matching family
membership and sequence inputs. Do not supply the baseline species tree in place
of each variant's inference. Also prepare both prespecified CPM variants, then
evaluate all six variants under the existing18-endpoint uncertainty protocol.

Machine-readable snapshot:ob_candidate_neighborhood_verified_20260916.json
SHA256:38a4cc8c16e8e45800c88af618f9b5a9b3d77ab6b3d0e5b79fa555b88e849514.
Twelve focused tests passed; full unit suite1,147 passed in27.64s.

## Downstream Native Admission

All four changed candidate arms completed the frozen inferred-tree pipeline and
passed independent native-output validation. All retain251,378 proteins, with
zero cross-candidate merges. The table reports group integrity, not accuracy.

| Arm | Job | Elapsed | Root HOGs | Supported constraints | Detached constraints |
| --- | --- | --- | ---: | ---: | ---: |
| norm_low | 21316_0 | 7:35 | 59,729 | 5,903 | 2,612 |
| norm_high | 21316_1 | 6:05 | 59,822 | 5,782 | 2,563 |
| margin_low | 21316_2 | 19:33 | 59,329 | 6,527 | 3,446 |
| margin_high | 21316_3 | 11:37 | 60,156 | 5,361 | 2,090 |

Per-arm ob_candidate_*_native_verified_20260916.json snapshots preserve scheduler,
command, native source/tool/input/tree/output and partition checks. Times are
Slurm task elapsed with input verification and reused checkpoints on a shared
machine, not end-to-end inference benchmarks. No parameter scores have been
calculated; the two CPM variants still require candidate/phylogeny execution
before the complete six-variant analysis.

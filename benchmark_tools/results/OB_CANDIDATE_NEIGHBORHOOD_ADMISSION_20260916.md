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

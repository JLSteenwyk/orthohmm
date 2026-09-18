# Sequence-Control Graph Memory Planning

## Finding

The frozen graph code memory-maps numeric checkpoint inputs but subsequently
allocates full-length boolean masks, selected hit copies, integer slot keys
and winner arrays. Therefore checkpoint loading is not evidence that the full
all-hit graph fits a given memory allocation. Do not silently truncate the
all-hit variant, and do not treat the independently prespecified top100
diagnostic as a replacement if all-hit execution fails.

Frozen `orthohmm/accuracy.py` SHA-256:
`1a35944ab7fea859143f599b1262aa111272787f6f735b9acc37d299427f2ad6`,
from core revision `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`.

## Allocation Calculation

`estimate_rbnh_array_payload.py` audits the checkpoint and verifies these
frozen source bytes. For finite native int32/int32/float64 hit arrays, let:

- H = total directed hits, E = nonself hits.
- S = genes times (maximum species ID plus one), matching native allocation.
- W = tied winning hit rows and K = occupied query/target-species slots.

At the source line immediately after allocation of `best_targets`, named
array payload is exactly `H + 24*E + 12*S + 8*W + 24*K` bytes. The tool
bounds this snapshot using `1 <= K <= W <= E` and `K <= S` when E is positive.
When no eligible hit exists, native code returns before this snapshot.
Input array logical size is reported separately as `16*H + 4*genes` bytes.
The model requires 64-bit NumPy indices. It does not measure or predict RSS.

Excluded costs include sorting temporaries, allocator overhead, gene strings,
input page residency, libraries, later edge deduplication, graph-library
construction, clustering, singleton assignment and refinement. In particular,
the snapshot upper bound is NOT an upper bound on total peak RAM. There is no
automatic graph-feasibility admission or recommended memory allocation.

## Completed OrthoBench Check

Both already-retained checkpoints were rehashed and numerically audited;
no graph inference was repeated. File identities were checked again after
the calculation. These observations are not extrapolated to QfO.

| Variant | Directed hits | Snapshot payload bounds, bytes | Input logical bytes, separate |
| --- | ---: | ---: | ---: |
| All hits | 100,099,147 | 2,532,653,619 to 3,403,835,787 | 1,602,591,864 |
| Post-search top100 | 48,991,663 | 1,254,966,519 to 1,717,288,815 | 784,872,120 |

Machine-readable reports: `ob_all_hits_rbnh_array_payload_20260918.json` and
`ob_top100_rbnh_array_payload_20260918.json`. Each binds its checkpoint,
numeric audit, individual input files, frozen core and estimator source.
27 focused tests pass, including eight randomized tied-score cases that
trace actual native named-array sizes at the specified line, early return,
large integer counts, malformed dimensions, wrong core bytes and a checkpoint
change after numeric audit.

## Corrected QfO Next Step

### Queued Admission-Gated Calculation

Job21798 is queued after successful numeric admission21792; the scheduler
confirms PENDING(Dependency). The frozen executor is
`928eba0430f2348647d30819d5a6fa23873cb19d` at
`benchmarks/work/publication_qfo_graph_payload_v1`. It requests2CPUs,64GiB,
4hours onbizon with no requeue. No graph inference is launched.

`plan_qfo_sequence_graph_memory.py` requires completed numeric-admission and
conversion accounting, exact frozen converter/admitter revisions, successful
source-equivalence status, both prespecified variant identities and the
recorded evidence file hashes. It runs the existing estimator against each
admitted checkpoint manifest hash and compares genes/hits/self-hit counts
with the numeric audit. It rechecks evidence before writing the report.
The estimator verifies the frozen native graph source and checkpoint contents.

Future output: `benchmarks/work/qfo_graph_payload_20260918.json`. An existing
output is refused. The report records both variant estimates and explicitly
sets graph feasibility, graph launch, accuracy evaluation and publication
readiness to false. A resource decision still requires review of actual
counts and costs beyond the modeled snapshot; top100 cannot replace all-hit
failure. No completed QfO payload calculation is claimed at submission.

42 focused tests pass, including8 new tests for both variants/exact hashes,
missing variant, path/count disagreements, false feasibility admission and
failed scheduler rejection before output creation. Batch shell syntax passes.

### Earlier Planning Observation

As of the recorded live check, DIAMOND job21789 remained running; its
incremental execution log recorded56 completed target searches out of78,
27,688,723,204 bytes of completed hit tables, and active zero-based target56.
These are progress observations, not a full-panel success claim or hit count.
Do not estimate the remaining workload by assuming equal-sized targets.

After independent numeric admission21792, apply the estimator separately to
the admitted all-hit and top100 checkpoints using their exact manifest hashes.
Review the actual numeric counts and file sizes together with available
resources and downstream graph costs before freezing the graph launcher.
The corrected HMM/DIAMOND hit-coverage diagnostic21793 remains separately
queued; it does not establish graph memory feasibility or matched sensitivity.
No graph job is submitted by this planning work, no running pipeline is
changed, and no new accuracy result is admitted.

Example after numeric admission, substituting an admitted checkpoint/hash
and a fresh output path:

```bash
python benchmark_tools/estimate_rbnh_array_payload.py --checkpoint CHECKPOINT --sha256 MANIFEST_SHA256 --core benchmarks/work/publication_method_native_v2/orthohmm/accuracy.py --output NEW_REPORT.json
```

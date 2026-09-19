# Corrected QfO Full Reconciliation: Native Validation

The final factorial configuration, `p1_c1_r1` (profile refinement,
candidate expansion and reconciliation enabled), completed under job
`21760_3` with exit `0:0` in `01:54:14`. Native admission job `21764`
completed with exit `0:0` in `00:02:23`.

The retained [admission receipt](qfo_corrected_factorial_native_admission_21764.json)
is 11,488 bytes, SHA-256
`c8a289ac8128711da6c4e93654ce6953e48854b8d21b2adfea074d3feb8c70aa`.
An additional recursive check independently rehashed all 23 distinct
file records directly embedded in the receipt, including the native pairs,
execution status, species tree, metrics and frozen helper sources. All
byte counts and SHA-256 values matched. This extra check is distinct from
the scheduled validator's full 145,649-artifact inventory check.

| Native integrity quantity | Value |
| --- | ---: |
| Candidate families | 351,739 |
| Candidate genes | 984,137 |
| Root HOGs | 366,068 |
| Root-HOG genes | 984,137 |
| Split source families | 7,234 |
| Cross-source merges | 0 |
| Native ortholog pairs | 5,959,560 |
| Membership constraints | 40,169 |
| Supported constraints | 28,835 |
| Detached constraints | 11,334 |
| Detached genes | 26,107 |
| Ortholog pairs removed by membership policy | 191,368 |
| Root HOGs added by membership policy | 11,186 |

All candidate genes are preserved. These are output-integrity and internal
pipeline counts, not truth-based accuracy or mechanistic validation.
The inferred species-tree SHA-256 is
`24198609aaaa7d9bd1c7f9a6ba7ee160e3450572be4c6869933681099afff26e`.
Native pairs are 495,749,756 bytes, SHA-256
`22916193404283c4257e682529b7012c95a7fd3b14a500d0f086d76a988b65f6`.

Reference conversion job `21772` completed with exit `0:0` in `00:03:31`.
All 5,959,560 native pairs mapped without loss. The
[conversion receipt](qfo_corrected_factorial_pairs_21772.json) is 32,353
bytes, SHA-256
`d5fea582e753aa6efb1421ba3ec0959d03ab31bc5b4d3d778e7370816c633e03`.
An independent recursive check verified all 106 distinct directly embedded
file records. The conversion's fresh native-admission check exactly
reproduced the retained admission receipt. Both mapped pair files have
91,726,400 bytes and SHA-256
`e75b14f69637d6c56f9a87debbd9b6915e30bb4bbf462a0be30197910f39293d`.

Assessment `21787` is running; independent scoring admission `21788`
remains downstream. No score is admitted by this report.
Reconciliation-on evaluation must
use native pair predictions, not Root-HOG clique pairs. The complete
factorial uncertainty job `21894` remains dependent on scoring admission.

Elapsed time above is provenance from the shared host, not a controlled
comparative runtime. No frozen configuration, endpoint or input changed.
The DGX panel remains undisturbed and publication readiness is not claimed.

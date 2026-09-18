# Corrected Proteinortho Native Output and Pair Conversion

Inference job 21708 completed with exit 0:0 in scheduler elapsed 1:18:48.
Its native log ends with group/graph output messages and `All finished.`
The guarded runner recorded successful execution after its input/runtime
checks. Independent admission 21715 completed with exit 0:0 in 1:01.

Native admission report: `qfo_corrected_proteinortho_native_20260918.json`,
SHA-256 `d04054e8d48463cdd3fe55756003b849e4ea0612cb2d996f2060ff108727775c`.
Rechecked the admission source and all 334 bound records after completion.
The graph covers 3,003 species-pair sections across 78 input proteomes and
984,137 accessions, with 4,695,385 valid relations. Its SHA-256 is
`1bf0128f230d80d247aed1e7a3f97024feaa31904ee4596f438f5785a8cc145d`.
The group table is preserved and hash-bound but is not the QfO pair input.

## Conversion

Conversion job 21717 completed with exit 0:0 in 34 seconds using detached
executor `publication_qfo_corrected_comparator_pairs_v1`, revision
`01104e032c0b5ffc704eedfc69dd820d96cd6d11`. The reviewed native-admission
checksum above was supplied explicitly. Shell syntax validation passed for
`qfo_corrected_comparator_pairs_batch_20260918.sh`.

| Quantity | Count |
| --- | ---: |
| Converted native pairs | 4,695,385 |
| Reference-mapped pairs | 4,695,385 |
| Mapping losses | 0 |
| Native duplicate relations | 0 |

Manifest: `qfo_corrected_proteinortho_pairs_20260918.json`, SHA-256
`c0620eb8772a983956d2f3b5cb98636bbaca2ba0c2adbcd9ee0ff0805d02342b`.
Raw and filtered files are both 71,089,234 bytes, SHA-256
`f37e39f0739c4900b63be473b8ef2fa33470c7712a5bff50aa5b6835da023a53`.
All conversion input/source/output records were rechecked after completion;
independent line counts equal 4,695,385 for both files. Large pair files
remain outside git. Participant ID is `qfo_corrected_proteinortho`.

## Resources and Limits

The retained GNU time record reports native wall time 1:18:29, user CPU
134,669.24 seconds, system CPU 1,923.05 seconds and maximum RSS 2,617,320
KiB. These are that wrapper's measurements, not a claim of aggregate
simultaneous process-tree memory. Shared-host inference is not controlled
matched-resource efficiency evidence and is separate from DGX timings.

Neither successful execution nor conversion establishes biological accuracy.
The complete native graph is the selected pair representation, not raw
search edges or inferred cliques from the group table. No inference was
repeated or original-release result overwritten.

## Progress Ledger

Previous turn: progress, seventh original factorial assessment retained as
`5f79868`. Current turn: reviewed corrected Proteinortho terminal execution
and admission, submitted/completed pinned conversion and verified zero-loss
mapping. Next: freeze and run all six QfO assessment endpoints for this
corrected participant, then independently admit outputs before reporting
scores. Other corrected tools, final factorial cell, uncertainty and
dedicated timing admission remain outstanding. Publication readiness has
not been established.

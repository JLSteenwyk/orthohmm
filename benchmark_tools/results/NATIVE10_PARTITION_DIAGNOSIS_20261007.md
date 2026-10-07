# Native10 Partition Diagnosis

Native inference job 23902 completed with exit 0:0. Its independent review
23910 failed with `Native partition loses or adds input genes`; the retained
failure has SHA256
`bad4de5092098cf461ad618a3b9527cf2f02800d35a09e37fe0e413447560e99`.
That failure remains authoritative for admission. Neither the diagnostic nor
the successful process exit replaces a terminal review.

After read-only preliminary comparisons, diagnostic source `c7ff87df` was
committed and pushed before one selected execution. Its 31 invented fixtures
passed. The source SHA256 is
`137062e9141fe5db7c07de6749ac5e4062d6bc0a046b253cfe6968911f264da2`.
The execution used Python 3.10.13, cleared Python/library injection variables,
disabled bytecode writing, fixed `PYTHONHASHSEED=0` and restricted numerical
libraries to one thread. The exact command is retained in the time record.

The report is `native10_partition_diagnosis_20261007_v1.json`, SHA256
`4d8ce7b56b07d9a14cc5ee511ca13ac7021cc865908b28fcd5db33c85f4bc4e1`.
It binds and rechecks 169 files, including all original/copied inputs, the
original failed receipt and unchanged parser/validator sources.

| Representation | Unique Genes | Memberships | Groups | Missing / Extra / Duplicate IDs |
|---|---:|---:|---:|---:|
| Checkpoint gene list | 984,137 | 984,137 | Not applicable | 0 / 0 / 0 |
| Final clustered file | 984,137 | 984,137 | 394,768 | 0 / 0 / 0 |
| Materialized orthogroups | 984,137 | 984,137 | 394,768 | 0 / 0 / 0 |
| Independently parsed root HOGs | 984,137 | 984,137 | 394,768 | 0 / 0 / 0 |
| Frozen-parser root HOGs | 984,137 | 984,137 | 394,768 | 0 / 0 / 0 |

All 78 per-species ownership/count checks match the retained preparation.
The checkpoint names are lexical. The independent and frozen parsers agree,
and the unchanged frozen root coverage gate passes. All four output partition
hashes are identical:
`52e4db26dedfb2a20c13b141ed7fb5a78b8e8eeb7f1f7bb5f89050c40c9ee70c`.

Reconstructing 392,110 source families from native labels gives exactly the
expected family IDs and canonical pre-phylogeny payload SHA256
`56624444cee3f79e92cd3da1bc054ba157f38fbdb33d19261b572c97877f358f`,
matching the native provenance manifest. This is a native-output consistency
check, not independent orthology truth.

The selected diagnostic exited 0 in 43.34 seconds, with 1,479,044 KiB peak RSS
and zero swaps. These are diagnostic costs on the shared host, not inference
timings or an efficiency claim. No alignment, search, clustering, tree inference,
resource replay, conversion or benchmark scoring was run.

Current retained coverage therefore contradicts the original review error.
Its cause remains **unestablished**. This does not prove a transient failure,
a parser defect, a pipeline defect, continuous output integrity or successful
full semantic/resource validation. All admission/next-identity flags remain
false. The next bounded investigation is the unchanged semantic kernel in the
original caller context, with a separate non-admitting diagnostic receipt.

## Original Caller Context

Context-probe source `49f5c35e` was committed/pushed before one selected
execution. Its source SHA256 is
`eb6db30bd691341f2189d3795a17b69b404c8a4f74087c760a3107ca0645f50e`.
The new receipt `native10_semantic_context_diagnosis_20261007_v1.json` has
SHA256 `d54f2062f1dd37a25a69c03b6e4993466f96ac2c5313692bc60f4580a732ab1c`.
It binds/rechecks 182 files, including the exact original request, amendment,
plan, baseline and caller/kernel sources. The original caller's context
construction is reproduced without changing any kernel function or globals.

The unchanged full semantic kernel passes: 984,137 input genes, 394,768
orthogroups/root HOGs and 5,115,410 native pair rows. Checkpoint array/hash,
frozen settings, materialization, source-family, species-tree and native-pair
checks all run through the existing kernel. This is stronger than the root
coverage-only probe, but still not the allocated terminal-review contract.
Its nested semantic success cannot authorize conversion, scoring or identity 11.
The original error's cause remains unestablished.

The context probe exited 0 in 58.10 seconds with 1,030,112 KiB peak RSS and
zero swaps. The exact command is in its time record. It did not rerun inference
or the completed resource/environment replay.

## Fresh Full Review Outcome

The prospectively documented, distinct full review23973 completed0:0 in44:35,
ended2026-10-07T14:27:43 in scheduler accounting. It invoked the unchanged
original reviewer with all runtime, resource, shared-environment and output
checks active. Its actual `review/review.json` SHA256 is
`7eab212d9788deb031232e4641089ff96503c0e5727e1ec369067f93ffb518b5`.
It binds55,162evidencefiles and has `status=native_success`,
`terminal_reviewed=true`, `native_outputs_validated=true`,
`primary_resources_replayed=true`, `shared_host_resources_reviewed=true` and
`next_identity_authorized=true`. Source is the original unchanged reviewer,
not a substituted diagnostic. These are terminal-review gates, not accuracy
assessment or publication readiness.

No inference was rerun. Primary retained native resources are53,995.367466727s
wall,1,585,693.069194CPU-s and19,355,951,104bytes step-lifetime peak memory,
under the precise original resource scopes. Shared-host distortion remains
unknown and potentially tool-dependent. Scientific timing/accuracy admission
remains false until the downstream contracts complete. The large7.43GB replay
and complete review stay outside Git. Original23910failure remains retained
with unestablished cause. The next authorized step is the frozen two-CPU QfO
pair conversion in a fresh namespace, then assessment and independent admission.

## Conversion Outcome

Distinct conversion23977 completed0:0 in2:44, ended14:33:56. Its frozen
converter receipt `native10_qfo_pairs_allocated_20261007_v1/results.json` has
SHA256 `593679a7c35a6c5e2eb706846a6fcced706dd0d624bccb5c1608e6f83f7e7040`.
It materialized5,115,410native inferred pairs; expected/total/retained counts
match and zero pairs were lost to the frozen QfO mapping. Both two-column files
are78,572,552bytes with SHA256
`6777d4c32d2f0cc22edc23006872c52be0d0ab6b8994b73b0f169be1b4c7fd90`.
All984,137inputs remain in the coverage denominator;542,336occur in at least
one predicted relation (fraction0.5510777462893885). This is relation coverage,
not accuracy, recall or complete orthology recovery. No benchmark endpoint has
yet been scored or independently admitted for this native attempt.

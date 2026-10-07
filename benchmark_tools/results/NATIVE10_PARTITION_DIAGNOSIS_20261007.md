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

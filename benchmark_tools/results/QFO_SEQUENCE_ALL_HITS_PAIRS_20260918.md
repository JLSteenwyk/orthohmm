# All-Hit Sequence Pair Conversion

Independent graph admission21823 completed before frozen conversion21824.
The conversion completed0:0 in2:34 with2CPUs/64GiB on the shared main host.
It used executor4f0c30e5cdf287a35c9600886aec0a41bcc0b720 and admission SHA-256
`ebd78f0602f1c57c1d3378b01f9179b9a83d677037d83ea0cc3035eafb766b30`.

The prespecified multipass_refined partition yielded11,300,151cross-species
group-derived clique pairs. Expected, emitted and retained counts agree;
zero pairs were removed by reference mapping. These are not native
phylogenetic ortholog predictions or search-graph edges.

- [Conversion report](qfo_sequence_all_hits_pairs_21824.json), SHA-256
  `afc9b35cdc9f648053015224abcf107d0373924292bac7fef82f14cb608ced33`.
- Pair payload:174,421,406bytes, SHA-256
  `19aa145724a231f1a20f7b22fec3ba440bb3e5f846d76646b0dfaf92aa588a76`.
  Retained under`benchmarks/results/qfo_sequence_pairs_v1/all_hits/`;
  large raw pairs are not committed.

Six-endpoint scoring21825 is RUNNING8CPUs/64GiB using frozen executor
425f0a7f5d5dc9e1438ab0a1766c45b13596dc14. Its preflight reports job21825,
variantall_hits andaccuracy_admitted=false. Queued independent admission21826
uses frozen executore798c159b9f65727d5b1055b98a81cba3a8f374d with2CPUs/64GiB,
afterany:21825. The validator requires successful terminal conversion and
scoring, reconstructs provenance and checks native outputs and six metrics.
An upstream failure cannot pass this gate; no automatic retry is used.

31focused assessment/admission/launcher tests pass. The new launcher rejects
missing/extra arguments, invalid variants and malformed job/revision IDs
before executor lookup. Existing tests exercise score provenance and output
validation. This is not a completed native assessment or evidence of matched
search sensitivity. No new accuracy score or controlled timing is claimed.

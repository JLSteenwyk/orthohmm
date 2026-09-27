# Canonical Candidate Control Result

Job 22323 completed `0:0` in 1:35. The
[independent readback](ob_canonical_candidate_readback_22323.json) validates
the pinned wrapper/driver plan, runtime identities, submission and terminal
scheduler result, reconstructed input-array identities, and all five candidate
partitions and merge/seed sidecars. Production inference remains unchanged.

The input auditor loads only the hash-verified admitted historical cache and
fresh numeric checkpoint. It independently sorts directed keys encoded as
`query * gene_count + target`, rather than calling the adapter's lexsort
policy, and reproduces all array dtype/shape/content hashes. The historical
and fresh score-vector hashes remain distinct: they were not rounded or
silently collapsed. Self-hit accounting remains 251,135; all four nonself
arms contain 18,235,373 directed keys.

All five canonical arms cover the same 251,378 genes exactly once and yield
**54,445 candidate families**. All ten pairwise partition comparisons are
identical. Every candidate file has SHA-256
`44f201dbc7e4eccf06d9401e5ce1b20fdf6ad523e15ee6007fe1c58f3941842d`,
the exact historical candidate-file hash.

| Canonical input arm | Genes changed versus historical candidates | Genes changed versus uncanonicalized fresh-order candidates |
|---|---:|---:|
| Historical order, historical scores | 0 | 93 |
| Historical order, fresh scores | 0 | 93 |
| Fresh order, historical scores | 0 | 93 |
| Fresh order, fresh scores | 0 | 93 |
| Full fresh self-hit control | 0 | 93 |

Thus canonical directed-pair ordering removes the observed order sensitivity
for these fixed seeds, score vectors, gene indexing and runtime. It happens
to restore the historical candidate file; the policy was defined and committed
before this full-data result and was not selected using benchmark accuracy.
This is not a cross-runtime or universal determinism result, nor a completed
production fix. Independent merge reconstruction establishes consistency with
each arm's candidate file, not equivalence of every merge constraint to the
historical trace or independent recomputation of search support.

Seventy-four focused tests pass, including independent sort/hash agreement,
invalid/duplicate-key rejection and a large unsigned-key regression preventing
implicit floating-point key conversion. Existing candidate-content tests are
included. The earlier startup failure, uncanonicalized controls and historical
scores remain retained. No full search, phylogeny, scoring or controlled timing
benchmark was performed by this experiment.

Next check the semantic membership constraints passed to reconciliation and
then validate downstream output under an explicitly identified experimental
configuration. Do not infer final F1 reproduction from candidate-file identity,
change historical score rows or promote a default without that evidence.

```bash
python -m benchmark_tools.audit_ob_canonical_candidates --repo . \
  --directory benchmarks/work/ob_canonical_candidates_20260926 \
  --output /tmp/canonical-candidate-readback.json
```

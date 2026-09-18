# Sampled Recorded-Event Pair Audit

For completed original-release factorial cell `p0_c1_r1`, independently
reconstructed pair predictions from the saved reconciliation-node table,
without calling OrthoHMM's pair-generation implementation. Selected the first
64 reconciled families ordered by SHA-256 of `20260923:family_id` (ID breaks
ties), from the original execution artifact inventory. No accuracy, reference
labels, group sizes or observed errors informed selection. This is a new
post-execution integrity audit, not a prospective accuracy endpoint.

The node/group tables must match the admitted execution's byte sizes and
SHA-256 values. Admission itself is pinned to
`a87eef1ae64a2d9eb8a39469761f08ee805e13ea8ae1332d1999ad72a3209ba9`.
The verifier checks rooted topology, complete descendant partitions, unique
leaf genes, descendant species unions and cross-child species overlap. Under
positive-paralogy semantics, overlap blocks cross-child ortholog pairs;
mapping conflict without overlap produces an uncertain event whose
cross-species pairs remain eligible. Expanded-cell membership filtering is
reconstructed by restricting those pairs to equal recorded final root-HOG IDs.

| Quantity | Count |
| --- | ---: |
| Eligible reconciled families | 24,320 |
| Deterministically sampled families | 64 |
| Sampled genes | 1,437 |
| Sampled nodes | 2,808 |
| Pairs reconstructed before membership filtering | 13,863 |
| Pairs excluded by recorded group boundaries | 725 |
| Retained reconstructed/native pairs | 13,138 |
| Missing or extra native pairs within sample | 0 |

Machine-readable report: `qfo_event_pairs_p0_c1_r1_20260918.json`, SHA-256
`02eddc676302f6cb2d5f41e848d7570181f13ffc0ff76632f0a4cc60ba22d4f6`.
Fifteen tests cover speciation, uncertain events, duplication with
cross-species child pairs, group filtering, malformed/disconnected trees,
incorrect labels, species sets and descendant coverage.

```sh
python benchmark_tools/audit_reconciled_pair_events.py \
  --admission benchmark_tools/results/qfo_factorial_native_p0_c1_r1_20260918.json \
  --admission-sha256 a87eef1ae64a2d9eb8a39469761f08ee805e13ea8ae1332d1999ad72a3209ba9 \
  --families 64 --output /tmp/qfo-recorded-event-pairs.json
```

## Scope Limits

This is not independent tree inference, species-tree reconciliation or
biological truth validation. Mapping-conflict annotations and final group
membership are conditioned on, not independently regenerated. Confidence
labels and which membership constraints should be accepted are not audited.
Bypass families and unsampled pairs are outside scope. Do not extrapolate
zero discrepancies in this deterministic sample to an error-rate confidence
interval or whole-output guarantee. Native scoring remains separately gated.

The sequence-search control review remains unchanged: equal E-value cutoffs,
hit caps, or matched hit counts are not proof of matched sensitivity. This
pair-conversion audit does not resolve that separate HMM-contribution question.

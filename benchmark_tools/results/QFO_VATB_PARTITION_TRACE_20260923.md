# Post-Hoc VATB Partition Trace

VATB was selected after observing the largest family-level F1 decline in the
low-CPM sensitivity result. This is an error-localization analysis, not a new
independent test or a justification for outcome-selected defaults.

The 28 represented reference accessions are identical in the control and low-CPM
raw SwissTrees counts. Native counts without priors change from TP=330, FN=27,
FP=0, TN=21 to TP=0, FN=357, FP=0, TN=21. Native scoring applies its family
pseudocount convention separately; these raw counts are not those adjusted
statistics.

## Retained Stage Membership

| Arm | Stage | Groups containing reference genes | Reference genes in singleton groups | Co-grouped reference pairs |
| --- | --- | ---: | ---: | ---: |
| Low CPM | Multipass | 2 | 0 | 351 |
| Low CPM | Multipass refined | 28 | 27 | 0 |
| Low CPM | Strict profiles | 2 | 0 | 351 |
| Low CPM | Strict profiles refined | 28 | 27 | 0 |
| Low CPM | Candidate expansion | 28 | 27 | 0 |
| Control | Strict profiles refined | 2 | 0 | 351 |
| Control | Candidate expansion | 2 | 0 | 351 |

Before each low-CPM refinement, the two groups contain 276 and 21 total genes;
27 reference genes occur together in the first group, and one in the second.
After refinement, the 27 reference genes from the first group are singletons.
The other reference gene remains in a 21-gene group. Candidate expansion does
not reunite any of these reference genes. The control's corresponding final
seed groups contain 189 and 21 genes, and its candidate groups contain 200 and
21 genes.

Co-grouped pair counts are all unordered within-reference membership pairs,
not predicted ortholog counts or true positives. They do not distinguish
orthology from paralogy and must not be scored as native pair predictions.

## Supported Interpretation

The loss of within-family co-grouping occurs at refinement, before phylogenetic
inference, and is already present at its candidate input. Thus phylogenetic
reconciliation is not the first point of separation for this family in the
low-CPM arm. This does not establish the exact conditional branch, explain why
candidate expansion fails to reunite members, or characterize other families.
The strict-profile and multipass snapshots are separate saved stages; their
sequence must not be interpreted as proof that HMMs recovered a previously
refined partition without reviewing the replay control flow.

The frozen refinement source includes copy-number-based splitting. Whether
its exact predicates fired here still requires species-count and call-path
evidence; larger group size alone is insufficient. No scientific code or
parameters were changed and no inference was rerun.

[Machine-readable memberships and input identities](qfo_vatb_partition_trace_20260923.json)
include all 28 accession-to-gene mappings and the seven selected partitions.
Their digests were checked against retained stage admissions or the frozen
factorial manifest. Full-partition integrity relies on those admissions; this
trace independently checks reference membership completeness and uniqueness.
Seven focused tests cover membership accounting, ambiguity and missing inputs.

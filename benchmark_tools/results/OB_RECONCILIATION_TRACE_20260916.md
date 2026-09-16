# Reconciliation And Constraint Trace

## Verified Scope

Reconstructed the250 candidate families touching all70 OrthoBench references:
114 have saved reconciliation nodes and136 follow the unambiguous bypass rule.
For each candidate, the reconstructed pre-constraint root groups and subsequent
satellite constraints exactly recover its native final partition. All193
constraints inside these candidates are included, even when their logged sides
do not themselves touch reference genes:121 have a native high-confidence
supporting pair and72 detach their source. This is not the genome-wide audit
of all8,440 constraints.

The reconstruction follows frozen `species_overlap` root rules and
`positive_paralogy` pair rules, including ancestor splits propagated from
descendant root duplications. It reads saved calls, not new gene-tree inference.
Recorded source, node-table, selected tree/checkpoint, species-tree, input and
output identities were checked;564 provenance records passed final validation.
Fifteen focused tests include direct agreement with the native reconciler on
binary and multifurcating examples. Full unit suite:1,107 tests passed.

Machine-readable evidence: `ob_reconciliation_trace_20260916.json`, SHA256
`4d5822db16b7c326ce2ac923ecd98395ea4bf9d0ef77cfa0a3cae4a9d1d0ebb9`.
The snapshot includes all candidate group memberships, triggering node IDs,
constraint evidence and every affected reference pair.

## Pair Dispositions

| Disposition | Within-Reference Pairs |
| --- | ---: |
| Already in different candidate families | 16,264 |
| Separated by pre-constraint root-lineage rules | 134 |
| Separated subsequently by satellite constraints | 1,441 |
| Retained in final root HOGs | 22,894 |
| Total | 40,733 |

This partitions the same40,733 pair rows as the verified stage trace. Of1,575
post-candidate losses,134 occur at root-lineage extraction and1,441 at the
subsequent membership-constraint step. Attribution follows execution order;
a pair already separated by root-lineage rules is not also counted as a
constraint loss. This is not an order-independent counterfactual decomposition.

| Reference Family | Root-Lineage Losses | Constraint Losses |
| --- | ---: | ---: |
| RefOG001 | 14 | 0 |
| RefOG008 | 12 | 0 |
| RefOG011 | 92 | 385 |
| RefOG013 | 0 | 10 |
| RefOG014 | 0 | 95 |
| RefOG021 | 0 | 806 |
| RefOG042 | 16 | 0 |
| RefOG054 | 0 | 5 |
| RefOG061 | 0 | 140 |
| Other61 reference families | 0 | 0 |

Four families have root-lineage losses, six have constraint losses, with
RefOG011 in both sets. Among the six prespecified illustrations, the95 lost
pairs inRefOG014 and806 inRefOG021 occur at the constraint step. RefOG005,
RefOG024, RefOG038 andRefOG067 have no post-candidate within-family losses.
These observations include neutral cases and do not select new examples.

## Interpretation Limits

These are raw within-reference co-membership counts, including within-species
pairs, low-certainty assignments and overlapping reference membership. They
are not official recall, independently validated orthology, or proof that
each lost pair should have been retained. The native label `high` is a rule
category, not calibrated biological confidence.

Most observed post-candidate reference-pair losses occur at the satellite
constraint step. This identifies an executed code mechanism, not a reason to
remove it: the earlier unconstrained OrthoBench control lowered precision,
and did not establish an F1 improvement. See
[unconstrained control](ORTHOBENCH_UNCONSTRAINED_RESULTS_20260916.md).

Saved rooting/mapping calls have not been established as correct evolutionary
history. Independent tree validation, duplication/fragment/domain annotation,
search rejection and profile-edge tracing, QfO confirmation and the biological
application remain open. No method default or publication claim was promoted.

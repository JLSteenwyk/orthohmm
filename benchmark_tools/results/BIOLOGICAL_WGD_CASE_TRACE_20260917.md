# Prespecified WGD Case Trace

All six prospectively chosen cases are retained. This is a descriptive trace,
not a new selection, parameter search or statistical test. The
[machine-readable trace](biological_wgd_case_trace_20260917.json) reconstructs
seven candidate families from saved reconciliation nodes, applies recorded
satellite constraints and requires exact agreement with the admitted final
root-HOG memberships. Five families have reconciliation nodes; two are
unambiguous families that bypass tree inference.

## Stage Results

Homolog coverage is the number of available non-S. cerevisiae pillar members
in the union of the two anchor groups. Counts below are out of six. A merged
group can have complete coverage without separating the paralogs.

| Pair | Candidate state | Candidate coverage | Preconstraint root coverage | Final coverage | Homologs outside final anchor groups |
| --- | --- | ---: | ---: | ---: | --- |
| YER132C / YGL197W | merged | 6 | 5 | 5 | Suva_7.66 |
| YDR122W / YLR096W | merged | 6 | 6 | 6 | none |
| YER059W / YIL050W | separated | 6 | 6 | 6 | none |
| YLR284C / YOR180C | merged | not evaluated | not evaluated | not evaluated | different reference pillars |
| YBR147W / YOL092W | merged | 6 | 5 | 5 | Suva_4.396 |
| YCL048W / YDR522C | merged | 6 | 3 | 3 | Skud_3.16, Smik_3.28, Suva_3.165 |

All six pairs have distinct anchor groups after root-lineage reconstruction,
before constraints, and in the final output. The reference-excluded pair is
not assigned a homolog-support or coverage score.

The five listed homologs are not missing from the output. Suva_7.66 is assigned
to RootHOG0005922 in Family0003871; Suva_4.396 to RootHOG0004005 in
Family0002367; and the three homologs of YCL048W/YDR522C to RootHOG0005062 in
Family0003151. All were present in the same candidate as the focal anchors.
Their separation from the anchors is reproduced before satellite constraints.
Thus candidate omission and later constraint splitting do not explain these
five observed coverage losses. This does not rule out earlier search effects
on candidate composition or resulting gene-tree topology.

Two selected families have unsupported logged constraints: event250 in
Family0000912 and event61 in Family0002367. Applying them does not change
the reported focal coverage or anchor-support outcomes. No claim is made
about all other families or the genome-wide effect of those constraints.

## Provenance and Limits

The candidate partition, merge log, reconciliation manifest and final root HOGs
are checked against the frozen native admission. Reviewed frozen phylogeny
source hashes and exact reconciliation rules are required. The species-tree
hash must match the manifest; each selected tree must match its retained
checkpoint, with the same candidate members and species-tree identity.
Every inspected artifact is rechecked after tracing.

The node table and family checkpoint files receive retrospective hashes in
this trace; they were not independently hash-admitted at native-run completion.
Exact final membership reconstruction and tree/checkpoint consistency provide
additional checks, not proof of an uninterrupted chain of custody.

The reconstruction uses saved node calls and the already tested frozen-rule
reconstruction helper. It does not reroot trees, independently infer the
evolutionary history, or distinguish an erroneous topology from an unsuitable
root-lineage rule. Reference pillars establish homology, not correct ancestral
copy assignments. No gene loss, functional difference or specific biological
mechanism is inferred from a group split alone.

## Reproduction

Requires retained native artifacts at their manifest paths; use a fresh output:

```sh
python benchmark_tools/trace_wgd_cases.py --repo . --output /tmp/wgd_case_trace.json
```

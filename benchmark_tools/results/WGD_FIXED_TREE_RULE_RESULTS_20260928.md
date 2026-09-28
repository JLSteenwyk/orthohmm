# Fixed-Tree Root-Rule Results

The [protocol](WGD_FIXED_TREE_RULE_PROTOCOL_20260928.md) was committed and pushed
as `28a1fc17` before alternative outcomes were evaluated. This is a post hoc
mechanism diagnostic on all six retained examples and seven candidate families,
not a new validation cohort or global accuracy comparison.

## Results

Entries are non-S. cerevisiae pillar homologs retained in either anchor group,
over the six available homologs. These are coverage counts, not F1 scores.

| Anchor pair | species_overlap | supported_children | confidence | mapped_event |
| --- | ---: | ---: | ---: | ---: |
| YER132C / YGL197W | 5/6 | 5/6 | 5/6 | 4/6 |
| YDR122W / YLR096W | 6/6 | 6/6 | 6/6 | 5/6 |
| YER059W / YIL050W | 6/6 | 6/6 | 6/6 | 6/6 |
| YLR284C / YOR180C | excluded | excluded | excluded | excluded |
| YBR147W / YOL092W | 5/6 | 5/6 | 5/6 | 4/6 |
| YCL048W / YDR522C | 3/6 | 3/6 | 3/6 | 1/6 |

The excluded case has anchors in conflicting reference pillars; it remains in
the output and is not silently replaced. All six anchor pairs remain separated
under every rule. Supported separation remains positive in four of the five
eligible cases, and negative for YCL048W / YDR522C, under every rule.

`supported_children` and `confidence` reproduce all seven baseline preconstraint
and final partitions exactly. `mapped_event` reduces focal coverage in four
eligible cases without improving supported separation. None of the alternatives
places any of the five previously traced homologs back with either anchor.
All genes remain assigned somewhere in the retained family partitions.

## Verification And Scope

The executable diagnostic first reproduces the retained trace and checks its
input provenance. The baseline matches all seven preconstraint/final partitions
and every biological scoring field for all six cases; only group display IDs
are disregarded. Ortholog pairs, paralog pairs and pair-confidence annotations
are identical across the four arms. Input hashes are rechecked after analysis.
Rooted trees, species tree, candidates, pair rule and logged constraints remain
fixed. Complete groups, case scores, five homolog destinations and artifact
hashes are retained in [the machine-readable result](wgd_fixed_tree_rules_20260928.json).

A fresh full replay matches every retained JSON field after JSON serialization
(in-memory group tuples serialize as lists). All 53 focused tests pass across
the diagnostic guards, trace reconstruction, independent score audit and WGD
scorer. These tests do not establish biological correctness of the trees.

These results rule out these three existing root-rule alternatives as a repair
for the five focal losses on these fixed inputs. They do not establish that
tree topology or rooting is erroneous, nor identify the true ancestral-copy
assignments. No method defaults, publication benchmark scores or uncertainty
intervals have changed. A wider intervention would need its own protocol and
validation; this selected-case analysis cannot justify a global recommendation.

Reproduce from the repository root with a fresh output path:

```bash
python -B -m benchmark_tools.probe_wgd_fixed_tree_rules --repo . --output /tmp/wgd_rules_replay.json
```
